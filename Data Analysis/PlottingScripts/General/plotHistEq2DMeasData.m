function plotHistEq2DMeasData
%plotHistEq2DMeasData Plot a histogram-equalized 2D data.

% Select a file to plot.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read the data file, convert the variable names, and specify the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

[~, base_filename] = fileparts(filename);

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev'))
        continue
    end
    dep_vals = data.(dep_name);
    dep_rels = data.rels.(dep_name);
    
    if isempty(dep_rels)
        disp(['Independent (sweep) variables for data variable ''',...
              strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end

    if length(dep_rels) == 1
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on only one sweep variable. ',...
              'This data will not be plotted.'])
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        if ~isempty(strfind(dep_name, 'Probability')) ||...
                ~isempty(strfind(dep_name, 'Amplitude')) ||...
                ~isempty(strfind(dep_name, 'Phase')) ||...
                strcmp(dep_name, 'I') ||...
                strcmp(dep_name, 'Q')
            
            dep_vals = (dep_vals - min(dep_vals(:))) ./...
                       (max(dep_vals(:)) - min(dep_vals(:)));
            dep_vals = histeq(dep_vals);
            extra_title = 'Line-by-Line Histogram-Equalized ';
            extra_filename = '_histeq_';
        else
            extra_title = '';
            extra_filename = '';
        end
        
        % Plot the data as a smooth surface.
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        zunits = getUnits(data, dep_name);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[extra_title, strrep(dep_name, '_', ' '), zunits, ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_smooth', extra_filename]));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[extra_title, strrep(dep_name, '_', ' '), zunits, ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_pixelated', extra_filename]));
    end
    if length(dep_rels) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
end

% Show a message box with the experiment parameters.
showMessageBox(data);