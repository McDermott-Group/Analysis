function data = loadMeasurementData(file)
%loadMeasurementData Load data from a data file.
%  DATA = loadMeasurementData(FILE) loads data from a data file FILE.
%  The function returns the loaded structure DATA.

if ~exist('file', 'var')
    % Select a single file to load.
    [filename, pathname, status] = selectMeasurementDataFile(1);
    if ~status
        data = struct();
        return
    end
    file = fullfile(pathname, filename);
end

% Read the data file, convert the variable names, and specify the units.
data = processMeasurementData(importMeasurementData(file));