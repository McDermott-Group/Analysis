function fitMeasData2Exponent(data_variable, data)
%fitMeasData2Exponent(DATA_VARIABLE, DATA) Fit data to an exponetial
%function, plot the data and the fit. 
%   data = fitMeasData2Exponent(DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to an exponent, and returns the data structure DATA with
%   the fit appended to it.

if ~exist('data_variable', 'var')
    error(['No dependent data variable to fit the exponent to is ',...
        'given as an input argument.'])
end

if ~exist('data', 'var')
    % Select a file.
    data = loadMeasurementData;
end
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

    f = fit(indep_vals(:), dep_vals(:),...
        'a * exp(-b * x) + c', 'StartPoint',...
        [dep_vals(1) - dep_vals(end),...
        1 / (max(indep_vals) - min(indep_vals) + 9 * eps),...
        dep_vals(end)]);
    
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
end