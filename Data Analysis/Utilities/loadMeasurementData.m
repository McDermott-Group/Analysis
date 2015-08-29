function data = loadMeasurementData
%loadMeasurementData Load data from a text data file into the MATLAB 
%workspace.

% Select a file to load.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read the data file, convert the variable names, and specify the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));