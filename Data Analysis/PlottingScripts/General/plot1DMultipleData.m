function plot1DMultipleData
%plot1DMultipleData Plot multiple 1D graphs in the same plot.

% Select files to plot.
[filenames, pathnames, status] = selectMeasurementDataFile;
if ~status
    return
end

if ~iscell(filenames)
    filenames = {filenames};
end
if ~iscell(pathnames)
    pathnames = {pathnames};
end

% Read the data file, convert the variable names, and specify the units.
try
    for k = 1:length(filenames)
        data{k} = processMeasurementData(importMeasurementData(fullfile(pathnames{k}, filenames{k})));
    end
catch
    error('The selected files contain unmatched data.')
end

% Create folder Plots in the same directory as the selected data file
% if it does not exist.
plts_path = fullfile(pathnames{1}, 'Plots');
if ~exist(plts_path, 'dir')
    mkdir(pathnames{1}, 'Plots')
end

for data_index = 1:length(data{1}.dep)
    dep_name = data{1}.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev'))
        continue
    end
    if isempty(data{1}.rels.(dep_name))
        disp(['Independent (sweep) variables for data variable ''',...
              strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end

    % Plot 1D data.
    if length(data{1}.rels.(dep_name)) == 1
        for k = 1:length(filenames)
            if ~isfield(data{k}, dep_name) || isempty(data{k}.rels.(dep_name))
                error('The selected files do not match.')
            end
        end
        indep_name = data{1}.rels.(dep_name){1};
        if length(data{k}.rels.(dep_name)) ~= 1 ||...
                ~strcmp(indep_name, data{k}.rels.(dep_name){1})
            error('The selected files do not match.')
        end
        if ~isfield(data{k}, indep_name)
            error('The selected files do not match.')
        end

        fields = fieldnames(data{1});
        legend_options = cell(length(fields), 1);
        for f = 1:length(fields)
            flag = true;
            if isnumeric(data{1}.(fields{f})) && length(data{1}.(fields{f})) == 1
                values = NaN(length(filenames), 1);
                for k = 1:length(filenames)
                    if ~isnumeric(data{k}.(fields{f})) || length(data{k}.(fields{f})) ~= 1
                        flag = false;
                        break   
                    else
                        values(k) = data{k}.(fields{f});
                    end
                end
            elseif ischar(data{1}.(fields{f}))
                values = cell(length(filenames), 1);
                for k = 1:length(filenames)
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
            if flag && length(filenames) == length(unique(values))
                legend_options{f} = fields{f};
            end
        end
        legend_options = legend_options(~cellfun('isempty', legend_options));
        if isempty(legend_options)
            error(['At least two files have indistinguishable sets of ',...
                'parameters. Is is likely that the same file was selected more than once.'])
        end
        legend_options_whitespaces = cellfun(@(str) strrep(str, '_', ' '),...
            legend_options, 'UniformOutput', false);
        choice = showLegendMenu(['Select a legend option for the ',...
            strrep(dep_name, '_', ' '), ' plot:'],...
            legend_options_whitespaces);
        if choice > 0
            legend_param = legend_options{choice};
            if isfield(data{1}.units, legend_param) && ~isempty(data{1}.units.(legend_param))
                lunits = [' ', data{1}.units.(legend_param)];
            else
                lunits = '';
            end
            legend_entries = cell(length(filenames), 1);
            for k = 1:length(filenames)
                if isnumeric(data{1}.(legend_param))
                    legend_entries{k} = [legend_options_whitespaces{choice}, ' = ',...
                        num2str(data{k}.(legend_param)), lunits];
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
        
        xmin = min(data{1}.(indep_name)(:));
        xmax = max(data{1}.(indep_name)(:));
        if length(data) > 1
            for k = 2:length(data)
                xmin = min([xmin, min(data{k}.(indep_name)(:))]);
                xmax = max([xmax, max(data{k}.(indep_name)(:))]);
            end
        end
        if xmax == xmin
            xmax = xmax + eps;
        end

        xunits = getUnits(data{1}, indep_name);
        yunits = getUnits(data{1}, dep_name);

        errorbar_flag = true;
        for k = 1:length(filenames)
            if ~isfield(data{k}, 'error') || ~isfield(data{k}.error, dep_name)
                errorbar_flag = false;
                break
            end
        end
        if errorbar_flag % Plot an errobar graph.
            createFigure('right');

            hold on
            for k = 1:length(filenames)
                errorbar(data{k}.(indep_name), data{k}.(dep_name),...
                    data{k}.error.(dep_name),...
                    '.', 'LineWidth', 1, 'MarkerSize', 15)
            end
            hold off
            xlim([xmin xmax])
            grid on

            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title(strrep(dep_name, '_', ' '), 'Interpreter', 'none', 'FontSize', 10)
            if choice > 0
                legend(legend_entries, 'Interpreter', 'none')
            end
            savePlot(fullfile(plts_path, [dep_name, '_errorbar']));
        end

        createFigure;
        hold on
        for k = 1:length(filenames)
            plot(data{k}.(indep_name), data{k}.(dep_name),...
                '.-', 'LineWidth', 1, 'MarkerSize', 15)
        end
        hold off
        xlim([xmin xmax])
        grid on

        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title(strrep(dep_name, '_', ' '), 'Interpreter', 'none', 'FontSize', 10)
        if choice > 0
            legend(legend_entries, 'Interpreter', 'none')
        end
        
        savePlot(fullfile(plts_path, [dep_name, '_simple']));
    end
    if length(data{1}.rels.(dep_name)) > 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than one sweep variable. ',...
              'The data will not be plotted.'])
    end
end