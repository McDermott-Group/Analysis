function plotDataVarVsDataVar(x_data_variable, y_data_variable)
%plotDataVarVsDataVar(X_DATA_VARIABLE, Y_DATA_VARIABLE)  Plot
%a scatter plot of a variable vs another one.
%   plotDataVarVsDataVar(X_DATA_VARIABLE, Y_DATA_VARIABLE) plots a scatter
%plot. X_DATA_VARIABLE is plotted along x axis and Y_DATA_VARIABLE is
%plotted along y axis.

if ~exist('x_data_variable', 'var') && ~exist('y_data_variable', 'var')
    error('Specify two dependent variable names as the input argument.')
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
[data, x_data_variable] = checkDataVar(data, x_data_variable);
[data, y_data_variable] = checkDataVar(data, y_data_variable);

x_vals = data.(x_data_variable)(:);
y_vals = data.(y_data_variable)(:);

if length(x_vals) ~= length(y_vals)
    error('The sizes are not equal for the specified data variables.')
end

xunits = getUnits(data, x_data_variable);
yunits = getUnits(data, y_data_variable);

createFigure;
plotDots(data.(x_data_variable), data.(y_data_variable))

xlabel([strrep(x_data_variable, '_', ' '), xunits], 'FontSize', 14);
ylabel([strrep(y_data_variable, '_', ' '), yunits], 'FontSize', 14);
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)

savePlot(fullfile(plts_path, [filename, '_',...
    x_data_variable, '-', y_data_variable]));