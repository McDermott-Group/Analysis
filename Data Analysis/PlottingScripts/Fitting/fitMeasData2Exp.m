function fitMeasData2Exp(data_variable, independent_variable)
%fitMeasData2Exp(DATA_VARIABLE, INDEPENDENT_VARIABLE) Fit data to 
%an exponetial function, plot the data and the fit. 
%DATA_VARIABLE should be a name of the data variable. INDEPENDENT_VARIABLE
% should be a name of the indepenedent variable.

% Select a file.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read the data file, convert the variable names, and specify the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

[~, base_filename] = fileparts(filename);

if ~exist('data_variable', 'var')
    error('Dependent data variable is not specified.')
end

data_variable = strrep(data_variable, ' ', '_');
if ~isfield(data, data_variable)
    error(['Data variable ''', strrep(data_variable, '_', ' '),...
        ''' is not found in the data.'])
end

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);

if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

if exist('independent_variable', 'var')
    found = false;
    for k = 1:length(dep_rels)
        if strcmp(dep_rels{k}, indepenedent_variable)
            found = true;
            break
        end
    end
    if ~found
        error('Dependent data variable is not specified.')
    end
% elseif length(dep_rels) > 1
%     error('More than one independent variable found.')
end

% Plot 1D data.
if length(dep_rels) == 1
    indep_name = dep_rels{1};
    indep_vals = data.(indep_name);

    f = fit(indep_vals(:), dep_vals(:),...
        'a * exp(-b * x) + c', 'StartPoint', [1, 0, 0]);
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    astr = ['a = ', num2str(f.a), ' ± ', num2str(ae), yunits];
%     be = max([abs(ci(1, 2) - f.b), abs(ci(2, 2) - f.b)]);
%     if ~isempty(xunits)
%         inv_xunits = [' (1/',  xunits(3:end)];
%     end
%     disp(['b = ', num2str(f.b), ' ± ', num2str(be), inv_xunits])
    be = max([abs(1 / ci(1, 2) - 1 / f.b), abs(1 / ci(2, 2) - 1 / f.b)]);
    invbstr = ['1/b = ', num2str(1 / f.b), ' ± ', num2str(be), xunits];
    ce = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    cstr = ['c = ', num2str(f.c), ' ± ', num2str(ce), yunits];
    
    if isfield(data, 'error') && isfield(data.error, data_variable) % Plot an errobar graph.
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
            plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ' [', data.Timestamp, ']'],...
               ['y(a,b,c,x) = a * exp(-b * x) + c: ', astr],...
               [invbstr, '; ' cstr]}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', data_variable, '_errorbar_expfit']));
    end

    createFigure;
    plotSimple(indep_vals, dep_vals, '.')  % Plot a simple 1D graph.
    hold on
        plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
    hold off
    legend('data', 'fit')
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14);
    title({[filename, ' [', data.Timestamp, ']'],...
           ['y(a,b,c,x) = a * exp(-b * x) + c: ', astr],...
           [invbstr, '; ' cstr]}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [base_filename, '_', data_variable, '_simple_expfit']));
end