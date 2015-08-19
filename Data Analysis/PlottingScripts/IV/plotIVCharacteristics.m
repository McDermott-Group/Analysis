function plotIVCharacteristics
%plotIVCharacteristics plots a current-voltage characteristics.

% Retrive the path that was used last time.
fid = fopen(fullfile(tempdir, 'plotIVCharacteristics_last_pathname.txt'), 'r');
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

[temp_filename, temp_pathname] = uigetfile({'*.ivc', 'Text Files';...
          '*.*', 'All Files'}, 'Select an IV characteristic',...
          fullfile(pathname, 'JPM.ivc'));

if isnumeric(temp_filename)
    return
end

% Save the current path.
fid = fopen(fullfile(tempdir, 'plotIVCharacteristics_last_pathname.txt'), 'w');
if fid ~= -1
    fprintf(fid, '%s', pathname);
    fclose(fid);
end

if ~exist([temp_pathname, temp_filename], 'file')
    disp(['Error in ', mfilename,': File ', temp_pathname, temp_filename, ' does not exist.']);
    return
end

data = importIVCharacteristic([temp_pathname, temp_filename]);
file = dir([temp_pathname, temp_filename]);
timestamp = file.date;

if isempty(data) || ~isfield(data, 'Current') || isempty(data.Current) ||...
        ~isfield(data, 'Voltage') || isempty(data.Voltage)
    error(['Error using ', mfilename('fullpath'), ': The selected file does not contain any meaningful data']);
end

createFigure;
plotSimple(1e3 * data.Voltage, 1e3 * data.Current);
xlabel('Voltage (mV)', 'FontSize', 14);
ylabel('Current (mA)', 'FontSize', 14);
title(['IV Characteristic [', timestamp, ']'], 'FontSize', 14)

% Saving the plot to Plots folder in the same directory as the selected
% data file
temp_plts_path = fullfile(temp_pathname, 'Plots');
if ~exist(temp_plts_path, 'dir')
    mkdir(temp_pathname, 'Plots')
end

[~, temp_filename] = fileparts(temp_filename);
savePlot(fullfile(temp_plts_path, temp_filename));