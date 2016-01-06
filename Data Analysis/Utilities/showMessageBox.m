function showMessageBox(data)
%showMessageBox Show meessage box with the experiment parameters.

fields = fieldnames(data);
params = cell(length(fields), 1);
for k = 1:length(fields)
    field = fields{k};
    value = data.(fields{k});
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
        params{k} = ['Comments: ', [value{:}]];
    elseif strcmp(field, 'Filename') || strcmp(field, 'Name') 
        continue
    elseif ischar(value)
        params{k} = [strrep(field,  '_', ' '), ': ', value];
    end
end
temp_struct.Interpreter = 'tex';
temp_struct.WindowStyle = 'non-modal';
msgbox(params(~cellfun('isempty', params)), 'Parameters', temp_struct)

end