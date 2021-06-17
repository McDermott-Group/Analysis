function plotMultiIQMeasData
%plotMultiIQMeasData   Plot multiple IQ space trajectories in the same plot.
%this is modified by LIU, for multiple runs

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

for data_index = 1:length(data{1}.dep)
    I_name = data{1}.dep{data_index};
    if ~contains(I_name, 'Is')
        continue
    end
    Q_name = strrep(I_name, 'I', 'Q');
    if ~isfield(data{1}, Q_name)
        continue
    end

    if length(data) > 1
        [legend_entries, choice] = selectLegendEntries(data,...
            [I_name, ' - ', Q_name]);
        if choice == 0
            continue
        end
        plot_title = {[strrep(I_name, '_', ' '), ' — ',...
            strrep(Q_name, '_', ' '), ' Trajectories'],...
            [strrep(filenames{1}, '_', '\_'), ' - ',...
            strrep(filenames{end}, '_', '\_')],...
            ['[', data{1}.Timestamp, ' - ',...
             data{end}.Timestamp, ']']};
    else
        [~, filename, ext] = fileparts(filenames{1});
        plot_title = [strrep(filename, '_', '\_'), ext,...
            ' [', data{1}.Timestamp, ']'];
    end

    xunits = getUnits(data{1}, I_name);
    yunits = getUnits(data{1}, Q_name);

    % Plot a simple trajectories.
    createFigure([.01, .1, .88, .8]);
    hold on
    
    for k = 1:length(data)
        scatter(data{k}.(I_name)(:), data{k}.(Q_name)(:), '.')
    end
    hold off
    grid on
    axis equal
    set(gca, 'box', 'on')

    xlabel([strrep(I_name, '_', ' '), xunits], 'FontSize', 14)
    ylabel([strrep(Q_name, '_', ' '), yunits], 'FontSize', 14)
    title(plot_title, 'FontSize', 10)
    if length(data) > 1
        legend(legend_entries, 'Interpreter', 'none', 'Location', 'Best')
    end

    savePlot(fullfile(plts_path, [I_name, '-', Q_name, '_simple']));
end