function plot2DMeanNormData(normalization_direction)
%plot2DMeanNormData(NORMALIZATION_DIRECTION) plots line-by-line mean normalized 2D data from a text
%data file. NORMALIZATION_DIRCTION should be either 'along_x' or 'along_y'.

if ~exist('normalization_direction', 'var') ||...
    (~strcmp(normalization_direction, 'along_x') &&...
     ~strcmp(normalization_direction, 'along_y'))
    normalization_direction = 'along_y'; 
end

% Select a file to plot.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read probability data file, convert the variable names, and define
% the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

% Create folder Plots in the same directory as the selected data file
% if it does not exist.
plts_path = fullfile(pathname, 'Plots');
if ~exist(plts_path, 'dir')
    mkdir(pathname, 'Plots')
end
[~, base_filename] = fileparts(filename);

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev'))
        continue
    end
    dep_vals = data.(dep_name);
    dep_rels = data.rels.(dep_name);
    
    if isempty(dep_rels) && print_messages
        disp(['Independent (sweep) variables for data variable ''', strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end

    if length(dep_rels) == 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on only one sweep variable. ',...
              'This data will not be plotted.'])
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        if ~isempty(strfind(dep_name, 'Probability')) ||...
                ~isempty(strfind(dep_name, 'ADC_Amplitude')) ||...
                strcmp(dep_name, 'I') || strcmp(dep_name, 'Q')
            if strcmp(normalization_direction, 'along_y')
                dep_vals = dep_vals ./ (mean(dep_vals, 2) * ones(1, size(dep_vals, 2)));
            elseif strcmp(normalization_direction, 'along_x')
                dep_vals = dep_vals ./ (ones(size(dep_vals, 1), 1) * mean(dep_vals));   
            end
            extra_title = 'Line-by-Line Max Normalized ';
            extra_filename = ['_mean_', normalization_direction];
        else
            extra_title = '';
            extra_filename = '';
        end
        
        % Plot the data as a smooth surface.
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[extra_title, strrep(dep_name, '_', ' '), ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(saveas(gca, fullfile(plts_path, [base_filename, '_', dep_name, '_smooth', extra_filename]));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[extra_title, strrep(dep_name, '_', ' '), ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(gca, fullfile(plts_path, [base_filename, '_', dep_name, '_pixelated', extra_filename]));
    end
    if length(dep_rels) > 2 && print_messages
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
end

% Show a message box with the experiment parameters.
showMessageBox(data);