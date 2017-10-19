function plotMeasData(data_variable, loadedData)
%plotMeasData   Plot data from a data file.
%   plotDataVar(DATA_VARIABLE) plots data from a selected data file.
%   If DATA_VARIABLE is specified, only the corresponding data will be
%   plotted. 

% Additional functionality Oct 2016: able to specify an already loaded data
%   structure that has been loaded via loadMeasurementData

% Select a file.
if ~exist('loadedData', 'var')
    data = loadMeasurementData;
else
    data = loadedData;
end

if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
    if isempty(dep_vars)
        return
    end
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
else
    plotDataVar(data, data_variable);
end

yeah=audioread('yeah.mp3');
sound(yeah,44e3);