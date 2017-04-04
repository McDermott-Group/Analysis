function fitMeasData2Lorentzian(data_variable)
%fitMeasData2Lorentzian(DATA_VARIABLE)  Fit data to a Lorentzian
%superimposed on a linear background.
%   fitMeasData2Lorentzian(DATA_VARIABLE) fits data for DATA_VARIABLE to
%   a Lorentzian superimposed on a linear background, plots the data and
%   the fit, and saves the data structure containing the fit.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable = selectDepDataVars(data, true);
    if isempty(data_variable)
        return
    end
    data_variable = data_variable{1};
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

if strcmp(data.units.(data_variable), 'dB')
    data.(data_variable) = 10.^(data.(data_variable) / 10);
    data.units.(data_variable) = 'Arb. Units';
    disp('The values in dB are converted to arbitrary power units.')
elseif strcmp(data.units.(data_variable), 'V')
    data.(data_variable) = data.(data_variable).^2;
    data.units.(data_variable) = 'V^2';
    disp(['The data values in V are squared to obtain a quantaty',...
         ' proportional to the power.'])
elseif strcmp(data.units.(data_variable), 'ADC Units')
    data.(data_variable) = data.(data_variable).^2;
    data.units.(data_variable) = 'ADC Units^2';
    disp(['The data values in ADC Units are squared to obtain a quantaty',...
         ' proportional to the power.'])
end

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);

if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

% Plot 1D data.
if length(dep_rels) == 1
    indep_name = dep_rels{1};
    indep_vals = data.(indep_name);

    f = BiasedLorentzian(indep_vals, dep_vals);
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    if strcmp(yunits, '')
        ampl_units = [xunits(1:end-1), '^2)'];
    else
        ampl_units = [yunits(1:end-1), '*', xunits(3:end-1), '^2)'];
    end
    amplitude_txt = ['Amplitude = ', num2str(f.a), ' ± ', num2str(ae),...
        ampl_units];
    be = max([abs(ci(1, 2) - f.b), abs(ci(2, 2) - f.b)]);
    f_c_txt = ['Resonance Frequency = ', num2str(f.b), ' ± ',...
        num2str(be),  xunits];
    ce = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    FWHM_txt = ['FWHM = ', num2str(abs(f.c)) ' ± ', num2str(ce), xunits];
    de = max([abs(ci(1, 4) - f.d), abs(ci(2, 4) - f.d)]);
    if strcmp(yunits, '')
        slp_units = [' (1/', xunits(3:end)];
    else
        slp_units = [yunits(1:end-1), '/', xunits(3:end)];
    end
    slope_txt = ['Slope = ', num2str(f.d),' ± ', num2str(de), slp_units];
    ee = max([abs(ci(1, 5) - f.e), abs(ci(2, 5) - f.e)]);
    background_txt = ['Background = ', num2str(f.e), ' ± ', num2str(ee),...
        yunits];

    Q = f.b ./ f.c;
    Q_txt = ['Q = ', num2str(Q), ' ± ', num2str(Q * (be / f.b + ce / f.c))];
    
    full_title = {[strrep(filename, '_', '\_'), ext,...
        ' [', data.Timestamp, ']'], f_c_txt, FWHM_txt, Q_txt, amplitude_txt,...
        [slope_txt, '; ', background_txt]};

    name = ['Fitted_', data_variable];
    data.units.(name) = data.units.(data_variable);
    data.rels.(name) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = name;
    data.plotting.(name).plot_title = full_title;
    data.(name) = f(indep_vals);

    % Plot an errobar graph.
    if isfield(data, 'error') && isfield(data.error, data_variable)
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
        plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(data_variable, '_', ' '), yunits], 'FontSize', 14)
        title(full_title, 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', data_variable,...
            '_lorentzianfit_errorbar']));
    end
    
    % Plot a simple 1D graph.
    createFigure;
    plotSimple(indep_vals, dep_vals, '.')
    hold on
    plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
    hold off
    legend('data', 'fit')
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
    ylabel([strrep(data_variable, '_', ' '), yunits], 'FontSize', 14)
    title(full_title, 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_lorentzianfit_simple']));

elseif length(dep_rels) == 2 % Plot 2D data.
    if ~isempty(strfind(dep_rels{1}, 'Frequency'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
    elseif ~isempty(strfind(dep_rels{2}, 'Frequency'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
        flip_fit = false;
    else
        error(['The data does not appear to depenend on any ',...
            '''Frequency'' variable.'])
    end

    dep_units = data.units.(data_variable);
    indep2_units = data.units.(indep_name2);

    f_c = zeros(length(indep_vals1), 3);
    FWHM = zeros(size(f_c));
    amplitudes = zeros(size(f_c));
    backgrounds = zeros(size(f_c));
    slopes = zeros(size(f_c));
    fitted = zeros(size(dep_vals));
    for k = 1:length(indep_vals1)
        f = BiasedLorentzian(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        amplitudes(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
        f_c(k, :) = [f.b, f.b - ci(2, 2), ci(1, 2) - f.b];
        FWHM(k, :) = [abs(f.c), f.c - ci(2, 3), ci(1, 3) - f.c];
        slopes(k, :) = [f.d, f.d - ci(2, 4), ci(1, 4) - f.d];
        backgrounds(k, :) = [f.e, f.e - ci(2, 5), ci(1, 5) - f.e];
        fitted(k, :) = f(indep_vals2);
    end
    if flip_fit
        fitted = fitted';
    end

    name = 'Extracted_Amplitude';
    data.(name) = amplitudes(:, 1);
    data.error.(name) = amplitudes(:, 2:3);
    if strcmp(dep_units, '')
        data.units.(name) = [indep2_units, '^2'];
    else
        data.units.(name) = [dep_units, '*', indep2_units, '^2'];
    end
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    
    name = 'Extracted_Resonance_Frequency';
    data.(name) = f_c(:, 1);
    data.error.(name) = f_c(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name)
    
    name = 'Extracted_FWHM';
    data.(name) = FWHM(:, 1);
    data.error.(name) = FWHM(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name)
    
    name = 'Extracted_Q';
    data.(name) = f_c(:, 1) ./ FWHM(:, 1);
    data.error.(name) = data.(name) .* ...
        (mean( f_c(:, 2:3), 2) ./  f_c(:, 1) +...
         mean(FWHM(:, 2:3), 2) ./ FWHM(:, 1));
    data.units.(name) = '';
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name)

    name = 'Extracted_Slope';
    data.(name) = slopes(:, 1);
    data.error.(name) = slopes(:, 2:3);
    if strcmp(dep_units, '')
        data.units.(name) = ['1/', indep2_units];
    else
        data.units.(name) = [dep_units, '/', indep2_units];
    end
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    
    name = 'Extracted_Background_Level';
    data.(name) = backgrounds(:, 1);
    data.error.(name) = backgrounds(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, data_variable, 'pixelated');

    name = ['Fitted_', data_variable];
    data.(name) = fitted;
    data.units.(name) = dep_units;
    data.rels.(name) = dep_rels;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'pixelated');
    
    saveMeasData(data, [filename, '_', data_variable, '_lorentzfit'])
end

end

function f = BiasedLorentzian(x, y)

median_y = median(y);
max_y = max(y);
min_y = min(y);
if max_y - median_y >= median_y - min_y
    background = min_y;
    idx = find(y == max_y, 1, 'first');
else
    background = max_y;
    idx = find(y == min_y, 1, 'first');
end

start = x(idx);
half_y = (median_y + y(idx)) / 2;
idx1 = find(y < half_y, 1, 'first');
idx2 = find(y < half_y, 1, 'last');
if ~isempty(idx1) && ~isempty(idx2) && x(idx1) ~= x(idx2)
    width = .1 * abs(x(idx2) - x(idx1));
else
    width = 0.01;
end

start_point = [.3 * (max_y - min_y) * width^2 / 4,...
    start, width, 0, background];

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.TolX = 1e-10;
opts.TolFun = 1e-10;
opts.StartPoint = start_point;

f = fit(x(:), y(:), '(a / ((x - b)^2 + (c / 2)^2)) + d * (x - b) + e', opts);

end