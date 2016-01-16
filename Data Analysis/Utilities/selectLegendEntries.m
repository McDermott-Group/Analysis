function [legend_entries, choice] = selectLegendEntries(data, dependent_variable)
%selectLegendEntries Select legend entry for multiple plots.
%   [LEGEND_ENTRIES, CHOICE] = selectLegendEntries(DATA, DEPENDENT_VARIABLE)
%   opens a legend selection diaglog (for varaiable DEPENDENT_VARIABLE
%   if the DEPENDENT_VARIABLE is specified). The function
%   analyzes the differences between the data structures inside
%   cell DATA. The output LEGEND_ENTRIES is a cell that contains legend
%   strings for each structure in cell DATA correspondingly. CHOICE
%   is a selected legend option (CHOICE is equal to zero if nothing is
%   selected).


if ~exist('data', 'var') && ~iscell(data)
    error('A cell containing data structures is expected.') 
end

cell_size = length(data);

fields = fieldnames(data{1});
legend_options = cell(length(fields), 1);
for f = 1:length(fields)
    flag = true;
    if isnumeric(data{1}.(fields{f})) && length(data{1}.(fields{f})) == 1
        values = NaN(cell_size, 1);
        for k = 1:cell_size
            if isfield(data{k}, fields{f})
                if ~isnumeric(data{k}.(fields{f})) ||...
                        length(data{k}.(fields{f})) ~= 1
                    flag = false;
                    break   
                else
                    values(k) = data{k}.(fields{f});
                end
            else
                flag = false;
            end
        end
    elseif ischar(data{1}.(fields{f}))
        values = cell(cell_size, 1);
        for k = 1:cell_size
            if ~ischar(data{k}.(fields{f}))
                flag = false;
                break   
            else
                values{k} = data{k}.(fields{f});
            end
        end
    else
        flag = false;
    end
    if flag && cell_size == length(unique(values))
        legend_options{f} = fields{f};
    end
end
legend_options = legend_options(~cellfun('isempty', legend_options));
if isempty(legend_options)
    error(['At least two files have indistinguishable sets of ',...
            'parameters. Is is likely that the same file was ',...
            'selected more than once.'])
end
legend_options_whitespaces = cellfun(@(str) strrep(str, '_', ' '),...
    legend_options, 'UniformOutput', false);
choice = showLegendMenu(['Select a legend option for the ',...
    strrep(dependent_variable, '_', ' '), ' plot:'],...
    legend_options_whitespaces);

legend_entries = cell(cell_size, 1);
if choice > 0
    legend_param = legend_options{choice};
    if isfield(data{1}.units, legend_param) &&...
            ~isempty(data{1}.units.(legend_param))
        lunits = [' ', data{1}.units.(legend_param)];
    else
        lunits = '';
    end
    for k = 1:cell_size
        if isnumeric(data{1}.(legend_param))
            legend_entries{k} = [legend_options_whitespaces{choice},...
                ' = ', num2str(data{k}.(legend_param)), lunits];
        elseif ischar(data{1}.(legend_param))
            if strcmp(legend_param, 'Filename')
                [~, name, ext] = fileparts(data{k}.(legend_param));
                legend_entries{k} = [name, ext];
            else
                legend_entries{k} = data{k}.(legend_param);
            end
        end
    end
end