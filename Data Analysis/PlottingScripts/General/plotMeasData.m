function plotMeasData(data_variable)
%plotMeasData   Plot data from a data file.
%   plotDataVar(DATA_VARIABLE) plots data from a selected data file.
%   If DATA_VARIABLE is specified, only the corresponding data will be
%   plotted.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
    for data_index = 1:length(dep_vars)
        dep_name = dep_vars{data_index};
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
else
    plotDataVar(data, data_variable);
end