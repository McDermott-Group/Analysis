function plotMulti1DMeasData(data_variable)
%plotMulti1DMeasData    Plot multiple 1D graphs in the same plot.
%   plotMulti1DMeasData(DATA_VARIABLE) plots DATA_VARIABLE from multiple
%   data sets in the same plot. Plots for all 1D data variable are
%   generated if DATA_VARIABLE is omitted.

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
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

for data_index = 1:length(data{1}.dep)
    dep_name = data{1}.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
            ~isempty(strfind(dep_name, '_Error')) ||...
            (exist('data_variable', 'var') &&...
            ~strcmp(dep_name, data_variable))
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
        indep_name = data{1}.rels.(dep_name){1};
        for k = 1:length(filenames)
            if  ~strcmp(indep_name, data{k}.rels.(dep_name){1}) ||...
                    ~isfield(data{k}, indep_name)
                error('The selected files do not match.')
            end
        end

        [legend_entries, choice] = selectLegendEntries(data, dep_name);
        if choice == 0
            continue
        end

        xunits = getUnits(data{1}, indep_name);
        yunits = getUnits(data{1}, dep_name);

        errorbar_flag = true;
        for k = 1:length(filenames)
            if ~isfield(data{k}, 'error') ||...
                    ~isfield(data{k}.error, dep_name)
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
            axis tight
            grid on

            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14);
            title(strrep(dep_name, '_', ' '), 'FontSize', 10)
            legend(legend_entries, 'Interpreter', 'none')
            savePlot(fullfile(plts_path, [dep_name, '_errorbar']));
        end

        % Plot a simple graph.
        createFigure;
        hold on
        for k = 1:length(filenames)
            plot(data{k}.(indep_name), data{k}.(dep_name),...
                '.-', 'LineWidth', 1, 'MarkerSize', 15)
        end
        hold off
        axis tight
        grid on

        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14);
        title(strrep(dep_name, '_', ' '), 'FontSize', 10)
        legend(legend_entries, 'Interpreter', 'none')
        
        savePlot(fullfile(plts_path, [dep_name, '_simple']));
    end
    if length(data{1}.rels.(dep_name)) > 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than one independent variable.'])
    end
end