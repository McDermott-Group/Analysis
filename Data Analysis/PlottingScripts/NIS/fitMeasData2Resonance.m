function fitMeasData2Resonance(data_variable)
%fitMeasData2Resonance(DATA_VARIABLE) Fit data to a resonance curve,
%plot the data and the fit. DATA_VARIABLE should be a name of the data
%variable.

if ~exist('data_variable', 'var')
    data_variable = 'Amplitude';
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

if strcmp(data.units.(data_variable), 'dB')
    data.(data_variable) = 10.^(data.(data_variable) / 20);
    data.units.(data_variable) = 'Arb. Units';
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

    [params, model] = skewfitquarterwave(indep_vals, dep_vals);
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    f0 = ['Resonance Frequency = ', num2str(params(1)), xunits, ';'];
    Qs = ['Internal Q = ', num2str(params(2)),...
        '; Coupling Q = ', num2str(params(3)), ';'];
    mismatch_A0 = ['Mismatch = ', num2str(params(4)),...
        '; Normalization Amplitude = ', num2str(params(5)), yunits];
    dep_vals = dep_vals / params(5);
    
    [~, fit] = model(params);
    if isfield(data, 'error') && isfield(data.error, data_variable)
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
            plot(indep_vals, fit, 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14)
        title({[filename, ext, ' [', data.Timestamp, ']'],...
            f0, Qs, mismatch_A0}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path,...
            [filename, '_', data_variable, '_resfit_errorbar']));
    end

    createFigure;
    plotSimple(indep_vals, dep_vals, '.')  % Plot a simple 1D graph.
    hold on
        plot(indep_vals, fit, 'r', 'Linewidth', 2)
    hold off
    legend('data', 'fit')
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
    ylabel('|S_{21}|', 'FontSize', 14)
    title({[filename, ext, ' [', data.Timestamp, ']'],...
        f0, Qs, mismatch_A0}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_resfit_simple']));
elseif length(dep_rels) == 2 % Plot 2D data.
    if strcmp(dep_rels{1}, 'RF_Frequency')
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
    elseif strcmp(dep_rels{2}, 'RF_Frequency')
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
    else
        error('The data does not appear to depenend on ''RF Frequency''.')
    end

    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    zunits = getUnits(data, data_variable);
    f0 = zeros(size(indep_vals1));
    Qi = zeros(size(f0));
    Qc = zeros(size(f0));
    L = zeros(size(f0));
    fit = zeros(size(dep_vals));
    for k = 1:length(indep_vals1)
        [params, model] = skewfitquarterwave(indep_vals2(:)', dep_vals(k, :));
        f0(k) = params(1);
        Qi(k) = params(2);
        Qc(k) = params(3);
        L(k) = params(4);
        dep_vals(k, :) = dep_vals(k, :) / params(5);
        [~, line_fit] = model(params);
        fit(k, :) = line_fit;
    end
    
    createFigure;
    plotSimple(indep_vals1, f0, '.-')
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
    ylabel(['Resonance Frequency ', yunits], 'FontSize', 14)
    title({[filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_resfreq_simple']));
    
    createFigure('right');
    plotSimple(indep_vals1, Qi, '.-')
    hold on
        plotSimple(indep_vals1, Qc, '.-')
    hold off
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
    ylabel('Quality Factor', 'FontSize', 14)
    legend('internal', 'coupling')
    title({[filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_qualfact_simple']));
    
    createFigure;
    plotSimple(indep_vals1, L, '.-')
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
    ylabel('Mismatch', 'FontSize', 14)
    title({[filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_rmismatch_simple']));
   
    createFigure;
    plotPixelated(indep_vals1, indep_vals2, fit')
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
    title({['|S21| Fit', zunits],...
        [filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_s21fit_pixelated']));
    
    createFigure('right');
    plotPixelated(indep_vals1, indep_vals2, dep_vals')
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
    title({['|S21| Data', zunits],...
        [filename, ext, ' [', data.Timestamp, ']']},...
        'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path,...
        [filename, '_', data_variable, '_s21data_pixelated']));
end
end

function [estimates, model] = skewfitquarterwave(f, y)
idx = find(y == min(y), 1, 'first');
start_f = f(idx);
half_max_y = (median(y) + y(idx)) / 2;
idx1 = find(y < half_max_y, 1, 'first');
idx2 = find(y < half_max_y, 1, 'last');
if ~isempty(idx1) && ~isempty(idx2) && f(idx1) ~= f(idx2)
    start_Qint = start_f / abs(f(idx2) - f(idx1));
else
    start_Qint = 20000;
end
model = @expfun;
start_point = [start_f, start_Qint, .8 * start_Qint, -.1, 1.1 * max(y)];
options = optimset('MaxFunEvals', 10000);
[estimates] = fminsearch(model, start_point, options);
    function [sse, signal] = expfun(params)
        f0 = params(1);
        Qi = params(2);
        Qc = params(3);
        L = params(4);
        A0 = params(5);
        Q = 1 / (1 / Qc + 1 / Qi);
        Smin = Qc / (Qi + Qc);
        dx = (f - f0) / f0;
        S21 = (Smin + 2 * 1i * Q * dx) ./ (1 + 2* 1i * Q * dx + 1i * L);
        signal = abs(S21);
        ErrorVector = signal - y / A0;
        sse = sum(ErrorVector.^2);
    end
end