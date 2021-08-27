function fitMeasData2Linear(data_variable, data)
%fitMeasData2Linear(DATA_VARIABLE, DATA)  Fit data to a linear
%
%   data = fitMeasData2Linear(DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to a linear function, and saves the data structure
%   containing the fit, plots the data and the fit, and saves the data
%   structure containing the fit.

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

    f = LinFit(indep_vals(:), dep_vals(:));
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    be = max([abs(1 / ci(1, 1) - 1 / f.p1), abs(1 / ci(2, 1) - 1 / f.p1)]);
    invbstr = ['1/b = ', num2str(-1 / f.p1, 4), ' ± ', num2str(be, 3),...
        xunits];
    ce = max([abs(ci(1, 2) - f.p2), abs(ci(2, 2) - f.p2)]);
    cstr = ['c = ', num2str(f.p2, 4), ' ± ', num2str(ce, 3), yunits];
    
    full_title = {[strrep(filename, '_', '\_'), ext,...
            ' [', data.Timestamp, ']'],...
            [strrep(data_variable, '_', ' '),' = - b * ',...
            strrep(indep_name, '_', ' '), ') + c'],...
            [astr, '; ', invbstr, '; ' cstr]};

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
            '_linearfit_errorbar']));
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
        '_linearfit_simple']));

elseif length(dep_rels) == 2 % Plot 2D data.
     % if ~isempty(strfind(dep_rels{1}, 'Delay')) ||...
     if ~isempty(strfind(dep_rels{1}, 'Time')) ||...
        ~isempty(strfind(dep_rels{1}, 'Duration')) ||...
        ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout')) ||...
        ~isempty(strfind(dep_rels{1}, 'QB_Drive_to_RO'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
     % elseif ~isempty(strfind(dep_rels{2}, 'Delay')) ||...
     elseif ~isempty(strfind(dep_rels{1}, 'Time')) ||...
             ~isempty(strfind(dep_rels{1}, 'Duration')) ||...
             ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout')) ||...
             ~isempty(strfind(dep_rels{1}, 'QB_Drive_to_RO'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
        flip_fit = false;
    else
        error(['The data does not appear to depenend on any ',...
            '''Time'', ''Duration'' or ''Qubit Drive to Readout'''...
            'variables.'])
    end

    dep_units = data.units.(data_variable);
    indep2_units = data.units.(indep_name2);

    time_const = zeros(length(indep_vals1), 3);
    rate_const = zeros(length(indep_vals1), 3);
    offset = zeros(size(time_const));
    fitted = zeros(size(dep_vals));
    for k = 1:length(indep_vals1)
        f = LinFit(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        time_const(k, :) = [-1/f.p1, 1 / (ci(1, 1)) - 1/f.p1,...
                                   1 / (f.p1) - 1/(ci(2, 1))];
        rate_const(k, :) = [-f.p1, (ci(1, 1) - f.p1),...
                                   (f.p1 - ci(2, 1))];
        offset(k, :) = [f.p2, f.p2 - ci(2, 2), ci(1, 2) - f.p2];
        fitted(k, :) = f(indep_vals2);
    end
    if flip_fit
        fitted = fitted';
    end
    
    name = 'Extracted_Time_Constant';
    data.(name) = time_const(:, 1);
    data.error.(name) = time_const(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    % plotDataVar(data, name, 'errorbar')
    
    name = 'Extracted_Rate_Constant';
    data.(name) = rate_const(:, 1);
    data.error.(name) = rate_const(:, 2:3);
    data.units.(name) = ['1/',indep2_units];
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')

    name = 'Extracted_Offset';
    data.(name) = offset(:, 1);
    data.error.(name) = offset(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
%     plotDataVar(data, name, 'errorbar')
    
    name = ['Fitted_', data_variable];
    data.(name) = fitted;
    data.units.(name) = dep_units;
    data.rels.(name) = dep_rels;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'pixelated');
    
    saveMeasData(data, [filename, '_', data_variable, '_linearfit'])
end
end

function f = LinFit(x, y)
%     opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts = fitoptions('Method', 'LinearLeastSquares',...
                      'Robust', 'LAR');
    f = fit(x(:), y(:), 'poly1', opts);
    
end