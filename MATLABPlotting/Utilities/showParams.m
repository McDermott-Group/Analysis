function showParams(data)
%showParams Show message box with the experiment parameters.

% Select a file.
if ~exist('data', 'var')
    data = loadMeasurementData;
    if isempty(fields(data))
        return
    end
end

flds = fieldnames(data);
params = cell(length(flds), 1);
for k = 1:length(flds)
    field = flds{k};
    value = data.(flds{k});
    if isnumeric(value)
        if length(value) == 1
            if isfield(data.units, field) && ~isempty(data.units.(field))
                units = [' ', data.units.(field)];
            else
                units = '';
            end
            params{k} = [strrep(field,  '_', ' '), ' = ',...
                num2str(value, 6), units];
        end
    elseif strcmp(field, 'Timestamp')
        params{k} = ['Time: ', data.(field)];
    elseif strcmp(field, 'Comments')
        params{k} = ['Comments: ', char(value{:})];
    elseif strcmp(field, 'Filename')
        continue
    elseif strcmp(field, 'Experiment_Name') ||  strcmp(field, 'Name') 
        params{k} = [strrep(field,  '_', ' '), ': ',...
            strrep(value, '_', '\_')];
   elseif strcmp(field, 'ExperimentName')
        params{k} = ['Experiment Name', ': ', strrep(value, '_', '\_')];
    elseif ischar(value)
        params{k} = [strrep(field,  '_', ' '), ': ', value];
    end
end
temp_struct.Interpreter = 'tex';
temp_struct.WindowStyle = 'non-modal';
neg=~cellfun(@isempty,params);
[~,ix] = sort(params(neg));   % sort the first column from 2nd row onwards and get the indices
temp=params(neg);
params(neg) = temp(ix);       % do not forget to add a 1 :)
msgbox(params(~cellfun('isempty', params)), 'Parameters', temp_struct)

end