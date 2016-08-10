function plotIVCharacteristics
%plotIVCharacteristics plots a current-voltage characteristics.

% Retrive the path that was used last time.
fid = fopen(fullfile(tempdir,...
    'plotIVCharacteristics_last_pathname.txt'), 'r');
if fid ~= -1
    pathname = fgetl(fid);
    fclose(fid);
end

% Open the user interface to select a file.
if ~exist('pathname', 'var')
    if exist('Z:\mcdermott-group\Data\Matched JPM', 'dir')
        pathname = 'Z:\mcdermott-group\Data\Matched JPM';
    elseif exist('Z:\Data\Matched JPM', 'dir')
        pathname = 'Z:\Data\Matched JPM';
    else
        pathname = '';
    end
end

[filename, pathname] = uigetfile({'*.ivc', 'Text Files';...
          '*.*', 'All Files'}, 'Select an IV characteristic',...
          fullfile(pathname, 'JPM.ivc'));

if isnumeric(filename)
    return
end

% Save the current path.
fid = fopen(fullfile(tempdir,...
    'plotIVCharacteristics_last_pathname.txt'), 'w');
if fid ~= -1
    fprintf(fid, '%s', pathname);
    fclose(fid);
end

if ~exist([pathname, filename], 'file')
    disp(['Error in ', mfilename, ': File ', pathname, filename,...
        ' does not exist.']);
    return
end

data = importIVCharacteristic([pathname, filename]);
file = dir([pathname, filename]);
timestamp = file.date;

if isempty(data) || ~isfield(data, 'Current') || isempty(data.Current) ||...
        ~isfield(data, 'Voltage') || isempty(data.Voltage)
    error(['Error using ', mfilename('fullpath'),...
        ': The selected file does not contain any meaningful data']);
end

% Plot raw data.
createFigure;
plotSimple(1e3 * data.Voltage, 1e3 * data.Current);
xlabel('Voltage (mV)', 'FontSize', 14)
ylabel('Current (mA)', 'FontSize', 14)
title(['Raw IV Characteristic [', filename, ', ', timestamp, ']'],...
       'Interpreter', 'none', 'FontSize', 10)

% Save the plot to Plots folder that is created in the same directory as
% the selected data file.
temp_plts_path = fullfile(pathname, 'Plots');
if ~exist(temp_plts_path, 'dir')
    mkdir(pathname, 'Plots')
end

[~, filename] = fileparts(filename);
savePlot(fullfile(temp_plts_path, [filename, '_raw']));

% Plot corrected data.
try
    data.Voltage = data.Voltage - mean(data.Voltage);

    [~, min_idx] = min(diff(data.Voltage));
    [~, max_idx] = max(diff(data.Voltage));

    min_I = data.Current(min_idx);
    max_I = data.Current(max_idx);

    V = data.Voltage(data.Current < max_I & data.Current > min_I);
    I = data.Current(data.Current < max_I & data.Current > min_I);

    cond = max(max(I)/max(V), min(I)/min(V));
    V_grnd = V(I ./ V > cond & abs(I) < .75 * min(max(I), abs(min(I))));
    I_grnd = I(I ./ V > cond & abs(I) < .75 * min(max(I), abs(min(I))));

    p = polyfit(I_grnd, V_grnd, 1);
    data.Voltage = data.Voltage - mean(V_grnd);

    R_app = (max(data.Voltage) - min(data.Voltage)) ./...
            (max(data.Current) - min(data.Current));

    createFigure;
    plotSimple(1e3 * (data.Voltage - p(1) * data.Current), 1e3 * data.Current);
    xlabel('Voltage (mV)', 'FontSize', 14)
    ylabel('Current (mA)', 'FontSize', 14)
    title({['Corrected IV Characteristic [', filename, ', ', timestamp, ']'],...
           ['Ground Resistance = ', num2str(p(1)), ' Ohm'],...
           ['Apparent Resistance = ', num2str(R_app), ' Ohm']},...
           'Interpreter', 'none', 'FontSize', 10)

    % Saving the plot to Plots folder in the same directory as the selected
    % data file
    temp_plts_path = fullfile(pathname, 'Plots');
    if ~exist(temp_plts_path, 'dir')
        mkdir(pathname, 'Plots')
    end

    [~, filename] = fileparts(filename);
    savePlot(fullfile(temp_plts_path, [filename, '_corrected']));
catch
    p = polyfit(data.Current, data.Voltage, 1);
    disp(['Apparent Resistance = ', num2str(p(1)), ' Ohm'])
end