function dep_vars = selectDepDataVars(data, single)
%selectDepDataVars Select a dependent data variable.
%   DEP_VARS = selectDepDataVars(DATA) opens a depenent data variable
%   selection diaglog. The output DATA_VARS is a cell that contains
%   selected dependent data variable names. Typically, this cell contains
%   either a single name or a list of all dependent of variables.


if ~exist('data', 'var')
    error(['A data structure or a cell containing data structures is '...
        'expected.']) 
end

if ~iscell(data)
    data = {data};
end
    
if ~exist('single', 'var')
    single = false;
end

cell_size = length(data);

fields = data{1}.dep;
dep_vars = cell(length(fields), 1);
for f = 1:length(dep_vars)
    include = true;
    for k = 1:cell_size
        found = false;
        for n = 1:length(data{k}.dep)
            if strcmp(fields{f}, data{k}.dep{n})
                found = true;
            end
        end
        if ~found
            include = false;
            break
        end
    end
    if include
        dep_vars{f} = fields{f};
    end
end
dep_vars = dep_vars(~cellfun('isempty', dep_vars));
if isempty(dep_vars)
    error('No (matching) dependent variables are found.')
end
dep_vars_menu = dep_vars;
if ~single
    dep_vars_menu{end+1} = 'All';
end
all_chosen = length(dep_vars) + 1;
dep_vars_menu = cellfun(@(str) strrep(str, '_', ' '),...
    dep_vars_menu, 'UniformOutput', false);
choice = showMenu('Select a dependent variable...',...
    dep_vars_menu);

if choice > 0
    if choice ~= all_chosen
        dep_vars = {dep_vars{choice}};
    end
else
    dep_vars = {};
end