function printOptimizationParams(index, data_variables)
%printOptimizationParams   Plot data from a data file.
%   printOptimizationParams(INDEX, DATA_VARIABLE) prints the optimization
%parameters for a given INDEX. If DATA_VARIABLES is specified, only
% the corresponding paremeter(s) will be printed.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variables', 'var')
    dep_vars = selectDepDataVars(data);
    if isempty(dep_vars)
        return
    end
else
    if ~iscell(data_variables)
        dep_vars = {data_variables};
    elseif iscell(data_variables)
        dep_vars = data_variables;
    else
        return
    end
end

if ~exist('index', 'var')
    if isfield(data, 'Contrast')
        [~, index] = max(data.Contrast);
    elseif isfield(data, 'Fidelity')
        [~, index] = max(data.Fidelity);
    elseif isfield(data, 'Occupation_Contrast')
        [~, index] = max(data.Occupation_Contrast);
    elseif isfield(data, 'Sideband_Suppression')
        [~, index] = max(data.Sideband_Suppression);
    end
end

for data_index = 1:length(dep_vars)
    dep_name = replace(dep_vars{data_index}, ' ', '_');
    depWname = replace(dep_name, '_', ' ');
    if contains(dep_name, '_Std_Dev') ||...
             contains(dep_name, '_Error')
        continue
    end
    vals = data.(dep_name);
    if index < length(vals)
        val = vals(index);
    end
    units = data.units.(dep_name);
    if strcmp(units, 'ns')
        val = round(val);
        format = '%d';
    elseif strcmp(units, 'dB')
        val = floor(4 * val) / 4;
        format = '%.2f';
    else
        format = '%f';
    end
    if ~isempty(units)
        prefix = ' * ';
    else
        prefix = '';
    end
    disp(['''', depWname, ''': ', num2str(val, format),...
        prefix, units, ',']) 
end

end