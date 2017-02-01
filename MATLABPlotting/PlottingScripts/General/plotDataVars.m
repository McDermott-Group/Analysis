function plotDataVars(data_variables)
%plotDataVars(DATA_VARIABLES)  Plot data variables listed in the cell
%in a sinle plot.
%   plotDataVars(DATA_VARIABLES) plots several data variables in a single
%   plot.

if ~iscell(data_variables)
    error('Data variable names should be fed in a cell.')
end

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);
plts_path = makeDirPlots(pathname);
plot_title = [filename, ext, ' [', data.Timestamp, ']'];

% Check that the data variable exists (compute it if necessary).
for k = 1:length(data_variables)
    [data, data_variable] = checkDataVar(data, data_variables{k});
    data_variables{k} = data_variable;
end

% Check that the data variables are all 1D arrays.
for k = 1:length(data_variables)
    if length(data.rels.(data_variables{k})) ~= 1
        error(['Data for variable ''', data_variables{k}, ''' is not 1D.'])
    end
end

indep = data.rels.(data_variables{1}){1};
xunits = getUnits(data, indep);
xlables = [strrep(indep, '_', ' '), xunits];
yunits = getUnits(data, data_variables{1});
ylables = data_variables;
legends = data_variables;

for k = 1:length(data_variables)
    if ~strcmp(indep, data.rels.(data_variables{k}){1})
        error(['The independent variables are not the same ',...
            'for specified data variables.'])
    end
    units = getUnits(data, data_variables{k});
    if ~strcmp(yunits, units)
        error('The dependent variable units do not match.')
    end
    ylables{k} = [strrep(data_variables{k}, '_', ' '), units];
    legends{k} = strrep(data_variables{k}, '_', ' ');
end

errorbar_flag = true;
for k = 1:length(data_variables)
    if ~isfield(data, 'error') ||...
            ~isfield(data.error, data_variables{k})
        errorbar_flag = false;
        break
    end
end

% Plot an errobar graph.
if errorbar_flag
    createFigure('right');
    hold on
    for k = 1:length(data_variables)
        indeps = data.rels.(data_variables{k});
        errorbar(data.(indeps{1}), data.(data_variables{k}),...
                 data.error.(data_variables{k}),...
            '.', 'LineWidth', 1, 'MarkerSize', 15)
    end
    hold off
    axis tight
    grid on
    set(gca, 'box', 'on')

    xlabel(xlables, 'FontSize', 14)
    ylabel(ylables, 'FontSize', 14)
    title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
    legend(legends, 'Interpreter', 'none', 'Location', 'Best')
    savePlot(fullfile(plts_path, [filename, '_comb_errorbar']));
end

% Plot a simple graph.
createFigure;
hold on
for k = 1:length(data_variables)
    if strcmp(data_variables{k},'Phase')
        data.(data_variables{k}) = unwrap(data.(data_variables{k}));
    end
        
    indeps = data.rels.(data_variables{k});
    plot(data.(indeps{1}), data.(data_variables{k}),...
        '.-', 'LineWidth', 1, 'MarkerSize', 15)
end
hold off
axis tight
grid on
set(gca, 'box', 'on')

xlabel(xlables, 'FontSize', 14)
ylabel(ylables, 'FontSize', 14)
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
legend(legends, 'Interpreter', 'none', 'Location', 'Best')
savePlot(fullfile(plts_path, [filename, '_comb_simple']));
end