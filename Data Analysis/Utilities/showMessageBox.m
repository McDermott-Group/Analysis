function showMessageBox(data)
%showMessageBox Show meessage box with the experiment parameters.

fields = fieldnames(data);
params = cell(length(fields), 1);
for k = 1:length(fields)
    if isnumeric(data.(fields{k}))
        if length(data.(fields{k})) == 1
            if isfield(data.units, fields{k}) && ~isempty(data.units.(fields{k}))
                units = [' ', data.units.(fields{k})];
            else
                units = '';
            end
            params{k} = [strrep(fields{k},  '_', ' '), ' = ',...
                        num2str(data.(fields{k})), ' ', units];
        end
    end
    if strcmp(fields{k}, 'Timestamp')
        params{k} = data.(fields{k});
    end
    if strcmp(fields{k}, 'Comments')
        params{k} = ['Comments: ', [data.(fields{k}){:}]];
    end
end
temp_struct.Interpreter = 'tex';
temp_struct.WindowStyle = 'non-modal';
msgbox(params(~cellfun('isempty', params)), 'Parameters', temp_struct)

end