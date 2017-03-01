function fitMeasData2Ramsey(data_variable, data)
%fitMeasData2Ramsey(DATA_VARIABLE, DATA)  Fit Ramsey data to a decaying
%sinusoid.
%   data = fitMeasData2Ramsey(DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to an exponent, plots the data and the fit, and saves
%   the data structure containing the fit.

if ~exist('data', 'var')
    % Select a file.
    data = loadMeasurementData;
end
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

    f = RamseyFit(indep_vals(:), dep_vals(:));
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    % ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    % astr = ['a = ', num2str(f.a, 4), ' ± ', num2str(ae, 3), yunits];
    
    % be = max([abs(ci(1, 2) - f.b), abs(ci(2, 2) - f.b)]);
    % bstr = ['b = ', num2str(f.b, 4), ' ± ', num2str(be, 3), yunits];
    
    T2e = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    T2str = ['T_2^* = ', num2str(f.c, 4), ' ± ', num2str(T2e, 3), xunits];
    
    TRe = max([abs(ci(1, 4) - f.d), abs(ci(2, 4) - f.d)]);
    TRstr = ['T_{Ramsey} = ', num2str(f.d, 4), ' ± ', num2str(TRe, 3),...
            xunits];

    % phe = max([abs(ci(1, 5) - f.e), abs(ci(2, 5) - f.e)]);
    % phstr = ['phase = ', num2str(f.e, 4), ' ± ', num2str(phe, 3)];

    full_title = {[strrep(filename, '_', '\_'), ext,...
            ' [', data.Timestamp, ']'],...
            [strrep(data_variable, '_', ' '),' = a + b * ',...'
            'e^{-(', strrep(indep_name, '_', ' ') ,'/ T_2^*)^2}',...
            'sin[2\pi(', strrep(indep_name, '_', ' '),...
            ' / T_{Ramsey}) + d]'],...
            [T2str, '; ', TRstr]};

    name = ['Fitted_', data_variable];
    data.units.(name) = data.units.(data_variable);
    data.rels.(name) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = name;
    data.plotting.(name).plot_title = full_title;
    data.(name) = f(indep_vals);
    
    % Plot an errorbar graph.
    if isfield(data, 'error') && isfield(data.error, data_variable)
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
        plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14)
        title(full_title, 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', data_variable,...
            '_expfit_errorbar']));
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
        '_expfit_simple']));

elseif length(dep_rels) == 2 % Plot 2D data.
    if ~isempty(strfind(dep_rels{1}, 'Duration')) || ...
            ~isempty(strfind(dep_rels{1}, 'Delay'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
    elseif ~isempty(strfind(dep_rels{2}, 'Duration')) || ...
            ~isempty(strfind(dep_rels{2}, 'Delay'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
        flip_fit = false;
    else
        error(['The data does not appear to depenend on any ',...
               '''Duration'' or ''Delay'' variables.'])
    end

    dep_units = data.units.(data_variable);
    indep2_units = data.units.(indep_name2);

    a = zeros(length(indep_vals1), 3);
    b = zeros(size(a));
    TRamsey = zeros(size(a));
    T2star = zeros(size(a));
    phase = zeros(size(a));
    fitted = zeros(size(a));
    for k = 1:length(indep_vals1)
        f = ExpFit(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        a(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
        b(k, :) = [f.b, f.b - ci(2, 2), ci(1, 2) - f.b];
        T2star(k, :) = [f.c, f.c - ci(2, 3), ci(1, 3) - f.c];
        TRamsey(k, :) = [f.d, f.d - ci(2, 4), ci(1, 4) - f.d];
        phase(k, :) = [f.e, f.e - ci(2, 5), ci(1, 5) - f.e];
        fitted(k, :) = f(indep_vals2);
    end
    if flip_fit
        fitted = fitted';
    end
    
    name = 'Extracted_a';
    data.(name) = a(:, 1);
    data.error.(name) = a(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    
    name = 'Extracted_b';
    data.(name) = b(:, 1);
    data.error.(name) = b(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    
    name = 'Extracted_T2_Star';
    data.(name) = T2star(:, 1);
    data.error.(name) = T2star(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    name = 'Extracted_Ramsey_Period';
    data.(name) = TRamsey(:, 1);
    data.error.(name) = TRamsey(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')

    name = 'Extracted_Phase';
    data.(name) = phase(:, 1);
    data.error.(name) = phase(:, 2:3);
    data.units.(name) = '';
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    
    name = ['Fitted_', data_variable];
    data.(name) = fitted;
    data.units.(name) = dep_units;
    data.rels.(name) = dep_rels;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'pixelated');
    
    saveMeasData(data, [filename, '_', data_variable, '_ramseyfit'])
end
end

function f = RamseyFit(x, y)
    ft = fittype('a + b * exp(-(x / c)^2) * sin(2 * pi * (x / d) + e)',...
        'independent', 'x', 'dependent', 'y' );
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    if y(1) < min(y(2:end)) || y(1) > max(y(2:end))
        x(1) = [];
        y(1) = [];
    end
    opts.Lower = [-max(abs(y)) -max(abs(y)) max(x) / 100 0 -2 * pi];
    opts.StartPoint = [y(1) max(y) - min(y) max(x) max(x) / 10 0];
    opts.Upper = [max(abs(y)) max(abs(y)) 10 * max(x) 10 * max(x) 2 * pi];
    f = fit(x, y, ft, opts);
end