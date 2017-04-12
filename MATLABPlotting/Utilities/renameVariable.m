function data = renameVariable(data, old_name, new_name)
%renameVariable Rename the data field.
%   DATA = renameVariable(DATA, OLD_NAME, NEW_NAME) renames OLD_NAME
%   variable to NEW_NAME in the data structure DATA and returns a new
%   data structure.
    if isfield(data, old_name)
        data.(new_name) = data.(old_name);
        data = rmfield(data, old_name);
    end
    if isfield(data, 'indep')
        for indep_index = 1:length(data.indep)
            if strcmp(data.indep{indep_index}, old_name)
                data.indep{indep_index} = new_name;
            end
        end
    end
    if isfield(data, 'dep')
        for dep_index = 1:length(data.dep)
            if strcmp(data.dep{dep_index}, old_name)
                data.dep{dep_index} = new_name;
            end
        end
    end
    if isfield(data, 'rels')
        rels_fields = fieldnames(data.rels); 
        for field_index = 1:length(rels_fields)
            if ~isempty(data.rels.(rels_fields{field_index}))
                indep_names = data.rels.(rels_fields{field_index});
                for indep_index = 1:length(indep_names)
                    if strcmp(indep_names{indep_index}, old_name)
                        indep_names{indep_index} = new_name;
                        data.rels.(rels_fields{field_index}) = indep_names;
                    end
                end
            end
        end
        if isfield(data.rels, old_name)
            data.rels.(new_name) = data.rels.(old_name);
            data.rels = rmfield(data.rels, old_name);
        end
    end
    if isfield(data, 'distr')
        if isfield(data.distr, old_name)
            data.distr.(new_name) = data.distr.(old_name);
            data.distr = rmfield(data.distr, old_name);
        end
    end
    if isfield(data, 'units')
        if isfield(data.units, old_name)
            data.units.(new_name) = data.units.(old_name);
            data.units = rmfield(data.units, old_name);
        end
    end
end