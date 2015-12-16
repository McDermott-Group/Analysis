function data = loadMeasurementData
%loadMeasurementData    Load data from a data file into the MATLAB 
%workspace. The function returns the loaded structure.

% Select a single file to load.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    data = struct();
    return
end

% Read the data file, convert the variable names, and specify the units.
file = fullfile(pathname, filename);
data = processMeasurementData(importMeasurementData(file));