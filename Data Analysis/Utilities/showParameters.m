function showParameters
%showParameters   how a message box with the experiment parameters.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

% Show a message box with the experiment parameters.
showMessageBox(data);