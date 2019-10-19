function fitMeasData2Exponent(data_variable, data)
%fitMeasData2Exponent(DATA_VARIABLE, DATA)  Fit data to an exponential
%function. 
%   data = fitMeasData2Exponent(DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to an exponent, and saves the data structure
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

    f = ExpFit(indep_vals(:), dep_vals(:));
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    astr = ['a = ', num2str(f.a, 4), ' ± ', num2str(ae, 3), yunits];
    be = max([abs(1 / ci(1, 2) - 1 / f.b), abs(1 / ci(2, 2) - 1 / f.b)]);
    invbstr = ['1/b = ', num2str(1 / f.b, 4), ' ± ', num2str(be, 3),...
        xunits];
    ce = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    cstr = ['c = ', num2str(f.c, 4), ' ± ', num2str(ce, 3), yunits];
    
    full_title = {[strrep(filename, '_', '\_'), ext,...
            ' [', data.Timestamp, ']'],...
            [strrep(data_variable, '_', ' '),' = a * exp(-b * ',...
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
    if ~isempty(strfind(dep_rels{1}, 'Time')) ||...
            ~isempty(strfind(dep_rels{1}, 'Duration')) ||...
            ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
    elseif ~isempty(strfind(dep_rels{2}, 'Time')) ||...
            ~isempty(strfind(dep_rels{2}, 'Duration')) ||...
            ~isempty(strfind(dep_rels{2}, 'Qubit_Drive_to_Readout'))
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
    amplitude = zeros(size(time_const));
    offset = zeros(size(time_const));
    fitted = zeros(size(dep_vals));
    for k = 1:length(indep_vals1)
        f = ExpFit(indep_vals2, dep_vals(k, :));
        ci = confint(f);
        amplitude(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
        time_const(k, :) = [1/f.b, 1 / ci(1, 2) - 1/f.b,...
                                   1 / f.b - 1/ci(2, 2)];
        offset(k, :) = [f.c, f.c - ci(2, 3), ci(1, 3) - f.c];
        fitted(k, :) = f(indep_vals2);
    end
    if flip_fit
        fitted = fitted';
    end
    
    name = 'Extracted_Amplitude';
    data.(name) = amplitude(:, 1);
    data.error.(name) = amplitude(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    name = 'Extracted_Time_Constant';
    data.(name) = time_const(:, 1);
    data.error.(name) = time_const(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')

    name = 'Extracted_Offset';
    data.(name) = offset(:, 1);
    data.error.(name) = offset(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    name = ['Fitted_', data_variable];
    data.(name) = fitted;
    data.units.(name) = dep_units;
    data.rels.(name) = dep_rels;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'pixelated');
    
    saveMeasData(data, [filename, '_', data_variable, '_expfit'])
end
end

function f = ExpFit(x, y)
%     opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts = fitoptions('Method', 'NonlinearLeastSquares',...
                      'Robust', 'LAR',...
                      'MaxFunEvals',10000,...
                      'MaxIter',10000,...
                      'DiffMaxChange',1,...
                      'Algorithm','Trust-Region');
    opts.Display = 'Off';
    opts.TolX = 1e-9;
    opts.TolFun = 1e-9;
    opts.StartPoint = [y(1) - y(end), 1 / (max(x) - min(x) + 9 * eps),y(end)];
    f = fit(x(:), y(:), 'a * exp(-b * x) + c', opts);
end