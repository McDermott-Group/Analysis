function plotMeasData
%plotMeasData   Plot data from a data file.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
             ~isempty(strfind(dep_name, '_Error'))
        continue
    end
    if ~isempty(strfind(dep_name, 'Phase'))
        data.(dep_name) = unwrap(data.(dep_name));
    end
    plotDataVar(data, dep_name);
end

% Show a message box with the experiment parameters.
showMessageBox(data);