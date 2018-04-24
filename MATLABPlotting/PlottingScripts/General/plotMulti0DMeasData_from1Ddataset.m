function plotMulti0DMeasData_from1Ddataset(data_variable)
%plotMulti0DMeasData_from1Ddataset    Plot multiple 0D graphs in the same
%plot from 1D datasets. For crunching the dummy variable

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
    data = cell(0);
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
elseif ~iscell(data_variable)
    dep_vars = {data_variable};
else
    dep_vars = data_variable;
end

dep_vars

for data_index = 1:length(dep_vars)
    dep_name = dep_vars{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
            ~isempty(strfind(dep_name, '_Error'))
        continue
    end
    if isempty(data{1}.rels.(dep_name))
        disp(['Independent (sweep) variables for data variable ''',...
              strrep(dep_name, '_', ' '), ''' are not specified.'])
    end

    % Plot 1D data.
    if length(data{1}.rels.(dep_name)) == 1
        for k = 1:length(filenames)
            if ~isfield(data{k}, dep_name) ||...
                    length(data{k}.rels.(dep_name)) ~= 1
                error('The selected files do not match.')
            end
        end
        indep_name = data{1}.rels.(dep_name){1}

        [legend_entries, choice] = selectLegendEntries(data, dep_name);
        if choice == 0
            continue
        end

        yunits = getUnits(data{1}, dep_name);

        errorbar_flag = true;
        for k = 1:length(filenames)
            if ~isfield(data{k}, 'error') ||...
                    ~isfield(data{k}.error, dep_name)
                errorbar_flag = false;
                break
            end
        end
        
        title_str_cell = {strrep(dep_name, '_', ' '),...
            [strrep(filenames{1}, '_', '\_'), ' - ',...
            strrep(filenames{end}, '_', '\_')],...
            ['[', data{1}.Timestamp, ' - ',...
             data{end}.Timestamp, ']']};


            hold off
            axis tight
            grid on
            set(gca, 'box', 'on')

            ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
            title(title_str_cell, 'FontSize', 10)
            legend(legend_entries, 'Interpreter', 'none', 'Location', 'Best')
            savePlot(fullfile(plts_path, [dep_name, '_errorbar']));
        end

        % Plot a simple graph.
        createFigure;
        hold on
        for k = 1:length(filenames)
            errorbar(k,mean(data{k}.Occupation),std(data{k}.Occupation),'.-', 'LineWidth', 1, 'MarkerSize', 15)
        end
        hold off
        axis tight
        grid on
        set(gca, 'box', 'on')


        ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
        title(title_str_cell, 'FontSize', 10)
        legend(legend_entries, 'Interpreter', 'none', 'Location', 'Best')
        
        savePlot(fullfile(plts_path, [dep_name, '_simple']));
    end
    if length(data{1}.rels.(dep_name)) > 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than one independent variable.'])
    end
end