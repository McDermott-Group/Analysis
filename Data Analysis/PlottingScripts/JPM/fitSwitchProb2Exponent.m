function data = fitSwitchProb2Exponent(data_variable, data)
%fitSwitchProb2Exponent(DATA_VARIABLE, DATA) Fit switching probability data
%to an exponential function, plot the data and the fit. 
%   data = fitMeasData2Exponent(DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to an exponent, and returns the data structure DATA with
%   the fit appended to it.

if ~exist('data', 'var')
    % Select a file.
    data = loadMeasurementData;
end
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable = selectDepDataVars(data, true);
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

    f = ExpFit(indep_vals(:), dep_vals(:));
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    astr = ['a = ', num2str(f.a, 4), ' ± ', num2str(ae, 3), ' 1/',...
        data.units.(indep_name)];
    
    full_title = {[strrep(filename, '_', '\_'), ext,...
            ' [', data.Timestamp, ']'],...
            [strrep(data_variable, '_', ' '),' = 1 - exp(-a * ',...
            strrep(indep_name, '_', ' '), ')'], astr};

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
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14);
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
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(data_variable, '_', ' '), yunits], 'FontSize', 14);
    title(full_title, 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_expfit_simple']));

elseif length(dep_rels) == 2 % Plot 2D data.
    if ~isempty(strfind(dep_rels{1}, 'Time'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
    elseif ~isempty(strfind(dep_rels{2}, 'Time'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
    else
        error(['The data does not appear to depenend on any ',...
            '''Time'' variable.'])
    end

    dep_units = data.units.(data_variable);
    indep2_units = data.units.(indep_name2);

    rate = zeros(length(indep_vals1), 3);
    for k = 1:length(indep_vals1)
        f = ExpFit(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        rate(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
    end
    
    name = 'Extracted_Rate';
    data.(name) = rate(:, 1);
    data.error.(name) = rate(:, 2:3);
    if strcmp(dep_units, '')
        data.units.(name) = ['1/', indep2_units];
    else
        data.units.(name) = [dep_units, '/', indep2_units];
    end
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    saveMeasData(data, [filename, '_', data_variable, '_expfit'])
end
end

function f = ExpFit(x, y)
f = fit(x(:), y(:),...
        '1 - exp(-a * x)', 'StartPoint', 1 / (max(x) - min(x) + 9 * eps));
end