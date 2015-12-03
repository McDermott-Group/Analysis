function plotMultiIQMeasData
%plotMultiIQMeasData   Plot multiple IQ space trajectories in the same plot.

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

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

for data_index = 1:length(data{1}.dep)
    I_name = data{1}.dep{data_index};
    if ~isempty(strfind(I_name, '_Std_Dev')) ||...
            isempty(strfind(I_name, 'I')) ||...
            length(data{1}.rels.(I_name)) ~= 1
        continue
    end
    Q_name = strrep(I_name, 'I', 'Q');
    if ~isfield(data{1}, Q_name)
        continue
    end

    for k = 1:length(data)
        if ~isfield(data{k}, I_name) ||...
                ~isfield(data{k}, Q_name) ||...
                length(data{k}.rels.(I_name)) ~= 1 ||...
                length(data{k}.rels.(Q_name)) ~= 1 ||...
                ~strcmp(data{1}.rels.(I_name){1}, data{k}.rels.(I_name){1}) ||...
                ~strcmp(data{1}.rels.(Q_name){1}, data{k}.rels.(Q_name){1}) ||...
                ~strcmp(data{k}.rels.(I_name){1}, data{k}.rels.(Q_name){1})
            error('The selected files do not match.')
        end
    end

    if length(data) > 1
        [legend_entries, choice] = selectLegendEntries(data,...
            [I_name, ' - ', Q_name]);
        if choice == 0
            continue
        end
        plot_title = [strrep(I_name, '_', ' '), ' — ',...
                      strrep(Q_name, '_', ' '), ' Trajectories'];
    else
        [~, filename, ext] = fileparts(filenames{1});
        plot_title = [filename, ext, ' [', data{1}.Timestamp, ']'];
    end

    xunits = getUnits(data{1}, I_name);
    yunits = getUnits(data{1}, Q_name);

    errorbar_flag = true;
    for k = 1:length(data)
        if ~isfield(data{k}, 'error') ||...
                ~isfield(data{k}.error, I_name) ||...
                ~isfield(data{k}.error, Q_name)
            errorbar_flag = false;
            break
        end
    end

    if errorbar_flag % Plot an errobar graph.
        createFigure([.9, .1, .88, .8]);
        marker_types = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
        hold on
        for k = 1:length(data)
            scatter(data{k}.(I_name), data{k}.(Q_name),...
                pi * 1.96^2 * data{k}.error.(I_name) .* data{k}.error.(Q_name),...
                linspace(1, 10, length(data{k}.(I_name))),...
                marker_types{mod(k - 1, length(marker_types)) + 1},...
                'filled', 'LineWidth', 1)
        end
        colormap(jet)
        hold off
        grid on
        axis equal

        xlabel([strrep(I_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(Q_name, '_', ' '), yunits], 'FontSize', 14);
        title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
        if length(data) > 1
            legend(legend_entries, 'Interpreter', 'none')
        end

        savePlot(fullfile(plts_path, [I_name, '-', Q_name, '_errorcirc']));
    end

    % Plot a simple trajectories.
    createFigure([.01, .1, .88, .8]);
    hold on
    for k = 1:length(data);
        plot(data{k}.(I_name), data{k}.(Q_name),...
            '.-', 'MarkerSize', 15)
    end
    hold off
    grid on
    axis equal

    xlabel([strrep(I_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(Q_name, '_', ' '), yunits], 'FontSize', 14);
    title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
    if length(data) > 1
        legend(legend_entries, 'Interpreter', 'none')
    end

    savePlot(fullfile(plts_path, [I_name, '-', Q_name, '_simple']));
end