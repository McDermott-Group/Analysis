function plot3DataVars(x_data_variable, y_data_variable, z_data_variable)
%plot3DataVars(X_DATA_VARIABLE, Y_DATA_VARIABLE, Z_DATA_VARIABLE)
%Plot a 3D trajectory.
%   plot3DataVars(X_DATA_VARIABLE, Y_DATA_VARIABLE, Z_DATA_VARIABLE) plots
%a 3D trajectory. X_DATA_VARIABLE is plotted along x axis, Y_DATA_VARIABLE
%is plotted along y axis, and Z_DATA_VARIABLE is plotted along z axis

if ~exist('x_data_variable', 'var') &&...
        ~exist('y_data_variable', 'var') &&...
        ~exist('z_data_variable', 'var')
    error('Specify three dependent variable names as the input argument.')
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
[data, z_data_variable] = checkDataVar(data, z_data_variable);

x_vals = data.(x_data_variable)(:);
y_vals = data.(y_data_variable)(:);
z_vals = data.(z_data_variable)(:);

len = length(x_vals);
if len ~= length(y_vals) || len ~= length(z_vals)
    error('The sizes are not equal for the specified data variables.')
end

xunits = getUnits(data, x_data_variable);
yunits = getUnits(data, y_data_variable);
zunits = getUnits(data, z_data_variable);

createFigure;
c = (1:len)';
surface([x_vals, x_vals], [y_vals, y_vals], [z_vals, z_vals], ...
    [c, c], 'EdgeColor','flat', 'FaceColor','none');
colormap(jet(len))
view(30, 45)
grid on

xlabel([strrep(x_data_variable, '_', ' '), xunits], 'FontSize', 14)
ylabel([strrep(y_data_variable, '_', ' '), yunits], 'FontSize', 14)
zlabel([strrep(z_data_variable, '_', ' '), zunits], 'FontSize', 14)
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)

savePlot(fullfile(plts_path, [filename, '_',...
    x_data_variable, '-', y_data_variable, '-', z_data_variable]));