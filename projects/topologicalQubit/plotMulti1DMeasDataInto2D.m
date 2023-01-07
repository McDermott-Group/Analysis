function plotMulti1DMeasDataInto2D(data_variable)
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
        
        title_str_cell = {[strrep(dep_name, '_', ' ') yunits],...
            [strrep(filenames{1}, '_', '\_'), ' - ',...
            strrep(filenames{end}, '_', '\_')],...
            ['[', data{1}.Timestamp, ' - ',...
            data{end}.Timestamp, ']']};
        
%         if errorbar_flag % Plot an errobar graph.
%             createFigure('right');
%             hold on
%             for k = 1:length(filenames)
%                 err = data{k}.error.(dep_name);
%                 if size(err, 1) == 2
%                     errorbar(data{k}.(indep_name), data{k}.(dep_name),...
%                         err(1, :), err(2, :),...
%                         '.', 'LineWidth', 1, 'MarkerSize', 15)
%                 elseif size(err, 2) == 2
%                     errorbar(data{k}.(indep_name), data{k}.(dep_name),...
%                         err(:, 1), err(:, 2),...
%                         '.', 'LineWidth', 1, 'MarkerSize', 15)
%                 else
%                     errorbar(data{k}.(indep_name), data{k}.(dep_name),...
%                         data{k}.error.(dep_name),...
%                         '.', 'LineWidth', 1, 'MarkerSize', 15)
%                 end
%             end
%             hold off
%             axis tight
%             grid on
%             set(gca, 'box', 'on')
%             
%             xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
%             ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
%             title(title_str_cell, 'FontSize', 10)
%             legend(legend_entries, 'Interpreter', 'none', 'Location', 'Best')
%             savePlot(fullfile(plts_path, [dep_name, '_errorbar']));
%         end
        
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
        set(gca, 'box', 'on')
        
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
        title(title_str_cell, 'FontSize', 10)
        legend(legend_entries, 'Interpreter', 'none', 'Location', 'Best')
        
        savePlot(fullfile(plts_path, [dep_name, '_simple']));
        
        %         createFigure;
        if length(filenames) > 1
            if choice >2
                for k = 1:length(filenames)
                    x = data{k}.(indep_name);
                    legend_str = cell2mat(legend_entries(k));
                    n = strfind(legend_str, '=');
                    y_dep = legend_str(1:n-1);
                    y(k) = sscanf(legend_str(n+1:end),'%f');
                    s = num2str(y(k));
                    ss = sscanf(s(end), '%d');
                    m = strfind(legend_str, num2str(ss));
                    yunits = legend_str(m+1:end);
                    z(k,:)=data{k}.(dep_name);
                    z1(k,:)=data{k}.(dep_name) - (ones(size(z(k,:), 1), 1) * median(z(k,:)));
                end
            else
                for k = 1:length(filenames)
                    x = data{k}.(indep_name);
                    y(k) = k;
                    legend_str = cell2mat(legend_entries(k));
                    y_dep = 'Measurement Sequence';
                    yunits = '#'
                    z(k,:)=data{k}.(dep_name);
                end
            end
            createFigure;
            imagesc(y, x, z')
            shading flat; view(0,90); colorbar; colormap jet; 
            set(gca,'YDir','normal');
            ylabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
            xlabel([y_dep, '(' yunits ')'], 'FontSize', 14)
            title(title_str_cell, 'FontSize', 10)
            
            savePlot(fullfile(plts_path, [dep_name, '_2D']));
            
            createFigure;
            imagesc(y, x, z1')
            shading flat; view(0,90); colorbar; colormap jet; 
            set(gca,'YDir','normal');
            ylabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
            xlabel([y_dep, '(' yunits ')'], 'FontSize', 14)
            title(['Cut Substracted ', title_str_cell], 'FontSize', 10)
            
            savePlot(fullfile(plts_path, [dep_name, '_Cut2D']));
        end
    end
    
    if length(data{1}.rels.(dep_name)) > 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
            ''' depends on more than one independent variable.'])
    end
end