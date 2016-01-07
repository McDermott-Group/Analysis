function fitMeasData2Lorentzian(data_variable)
%fitMeasData2Lorentzian(DATA_VARIABLE) Fit data to a Lorentzian curve,
%plot the data and the fit. DATA_VARIABLE should be a name of the data
%variable.

if ~exist('data_variable', 'var')
    error(['No dependent data variable to fit the Lorentzian to is ',...
        'given as an input argument.'])
end

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

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
    ampl = ['Amplitude = ', num2str(f.a), ' ± ', num2str(ae),...
        ampl_units, ';'];
    be = max([abs(ci(1, 2) - f.b), abs(ci(2, 2) - f.b)]);
    f0 = ['Resonance Frequency = ', num2str(f.b), ' ± ',...
        num2str(be),  xunits];
    ce = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    FWHM = ['FWHM = ', num2str(f.c) ' ± ', num2str(ce), xunits, ';'];
    de = max([abs(ci(1, 4) - f.d), abs(ci(2, 4) - f.d)]);
    if strcmp(yunits, '')
        slp_units = [' (1/', xunits(3:end)];
    else
        slp_units = [yunits(1:end-1), '/', xunits(3:end)];
    end
    slope = ['Slope = ', num2str(f.d),' ± ', num2str(de), slp_units, ';'];
    ee = max([abs(ci(1, 5) - f.e), abs(ci(2, 5) - f.e)]);
    background = ['Background = ', num2str(f.e), ' ± ', num2str(ee), yunits];

    % Plot an errobar graph.
    if isfield(data, 'error') && isfield(data.error, data_variable)
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
            plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ext, ' [', data.Timestamp, ']'],...
               f0, FWHM, ampl, slope, background},...
               'Interpreter', 'none', 'FontSize', 10)
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
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel(strrep(data_variable, '_', ' '), 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']'],...
           f0, FWHM, ampl, slope, background},...
           'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_lorentzianfit_errorbar']));
elseif length(dep_rels) == 2 % Plot 2D data.
    if ~isempty(strfind(dep_rels{1}, 'Frequency'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
    elseif ~isempty(strfind(dep_rels{2}, 'Frequency'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
    else
        error(['The data does not appear to depenend on any ',...
            '''Frequency'' variable.'])
    end

    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    dep_units = getUnits(data, data_variable);

    f0 = zeros(length(indep_vals1), 3);
    FWHM = zeros(size(f0));
    amplitudes = zeros(size(f0));
    backgrounds = zeros(size(f0));
    slopes = zeros(size(dep_vals));
    fit = zeros(size(dep_vals));
    for k = 1:length(indep_vals1)
        f = BiasedLorentzian(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        amplitudes(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
        f0(k, :) = [f.b, f.b - ci(2, 2), ci(1, 2) - f.b];
        FWHM(k, :) = [f.c, f.c - ci(2, 3), ci(1, 3) - f.c];
        amplitudes(k, :) = [f.d, f.d - ci(2, 4), ci(1, 4) - f.d];
        backgrounds(k, :) = [f.e, f.e - ci(2, 5), ci(1, 5) - f.e];
        fit(k, :) = f(indep_vals2);
    end
    
    createFigure;
    plotAssymErrorbar(indep_vals1, f0(:, 1), f0(:, 2), f0(:, 3))
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel(['Resonance Frequency ', yunits], 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']},...
            'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_resfreq_errbar']));
    
    createFigure('right');
    plotAssymErrorbar(indep_vals1, FWHM(:, 1), FWHM(:, 2), FWHM(:, 3))
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel(['FWHM ', yunits], 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']},...
            'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_fwhm_errbar']));
    
    createFigure('right');
    plotAssymErrorbar(indep_vals1, backgrounds(:, 1), backgrounds(:, 2),...
        backgrounds(:, 3))
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel(['Background ', dep_units], 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']},...
            'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_background_errbar']));
    
    createFigure('right');
    plotAssymErrorbar(indep_vals1, amplitudes(:, 1), amplitudes(:, 2),...
        amplitudes(:, 3))
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    if strcmp(dep_units, '')
        ampl_units = [yunits(1:end-1), '^2)'];
    else
        ampl_units = [dep_units(1:end-1), '*', yunits(3:end-1), '^2)'];
    end
    ylabel(['Amplitude ', ampl_units], 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']},...
            'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_amplitude_errbar']));
    
    createFigure('right');
    plotAssymErrorbar(indep_vals1, slopes(:, 1), slopes(:, 2), slopes(:, 3))
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    if strcmp(dep_units, '')
        slope_units = [' (1/', yunits(3:end)];
    else
        slope_units = [dep_units(1:end-1), '/', yunits(3:end)];
    end
    ylabel(['Slope', slope_units], 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']},...
            'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_errbar_errbar']));
   
    plotDataVar(data, data_variable);
    
    createFigure;
    plotPixelated(indep_vals1, indep_vals2, fit')
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
    title({['Lorentzian Fit with Linear Background to ',...
        strrep(data_variable, '_', ' '), dep_units],...
        [filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_lorentzianfit_pixelated']));
end
end

function f = BiasedLorentzian(x, y)

median_y = median(y);
max_y = max(y);
min_y = min(y);
if max_y - median_y >= min_y - median_y
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
    width = abs(x(idx2) - x(idx1));
else
    width = 0.01;
end

start_point = [.5 * (max_y - min_y) * width^2 / 4,...
    start, width, 0, background];
f = fit(x(:), y(:), '(a / ((x - b)^2 + (c / 2)^2)) + d * x + e',...
        'StartPoint', start_point);
end