function [data, data_variable] = checkDataVar(data, data_variable)
%checkDataVar   Check that the data variable exists in the data structure.
%DATA should be a data structure. DATA_VARIABLE should be a name of
%the data variable. The function returns structure DATA that
%contains the requested variable. This format is useful when the varible
%does not present in the original data structure but it is possible to
%compute it. The returned DATA_VARIABLE has its witespaces replaced by
%underscores.

if ~ischar(data_variable)
    error('The data variable name should be specified as a string.')
end
data_variable = strrep(data_variable, ' ', '_');
if ~isfield(data, data_variable)
    if ~isempty(strfind(data_variable, 'Phase_Space'))
        data = maximizeIQContrast(data);
    end
    if ~isfield(data, data_variable)
        error(['Data variable ''', strrep(data_variable, '_', ' '),...
            ''' is not found in the data.'])
    end
end

if isfield(data, 'dep')
    found = false;
    for k = 1:length(data.dep)
        if strcmp(data.dep{k}, data_variable)
            found = true;
            break
        end
    end
    if ~found
        error(['Data variable ''', strrep(data_variable, '_', ' '),...
        ''' is not found among the dependent data variables.'])
    end
end

end

