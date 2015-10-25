function plotMeasDataDiff
%plotMeasDataDiff Plot the difference between two text data sets.

% Select files to compute the difference.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file to subtract from...',...
     'Select the second data file to be subtracted...'});
if ~status
    return
end

% Read the data files, convert the variable names, and specify the units.
data1 = processMeasurementData(importMeasurementData(fullfile(pathnames{1}, filenames{1})));
data2 = processMeasurementData(importMeasurementData(fullfile(pathnames{2}, filenames{2})));

% Create folder Plots in the same directory as the selected data file
% if it does not exist.
plts_path = fullfile(pathnames{1}, 'Plots');
if ~exist(plts_path, 'dir')
    mkdir(pathnames{1}, 'Plots')
end
[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

for data_index = 1:length(data1.dep)
    dep_name = data1.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev'))
        continue
    end
    dep_vals1 = data1.(dep_name);
    if ~isfield(data2, dep_name)
        error('The selected files do not match.')
    else
        dep_vals2 = data2.(dep_name);
    end
    dep_rels1 = data1.rels.(dep_name);
    dep_rels2 = data2.rels.(dep_name);
     
    if (isempty(dep_rels1) || isempty (dep_rels2))
        disp(['Independent (sweep) variables for data variable ''',...
              strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end
    
    % Plot 1D data.
    if length(dep_rels1) == 1
        indep_name = dep_rels1{1};
        if length(dep_rels2) ~= 1 || ~strcmp(indep_name, dep_rels2{1})
            error('The selected files do not match.')
        end
        indep_vals = data1.(indep_name);
        if ~isfield(data2, indep_name) || length(indep_vals) ~= length(data2.(indep_name)) ||...
                any(indep_vals ~= data2.(indep_name))
            error('The selected files do not match.')
        end

        xunits = getUnits(data1, indep_name);
        yunits = getUnits(data1, dep_name);
        
        if isfield(data1, 'error') && isfield(data1.error, dep_name) &&...
           isfield(data2, 'error') && isfield(data2.error, dep_name) % Plot an errobar graph.
            difference = dep_vals1 - dep_vals2;
            difference_error = 1.96 * sqrt(data1.error.(dep_name).^2 + data2.error.(dep_name).^2);
            createFigure('right');
            plotErrorbar(indep_vals, difference, difference_error);
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title({[strrep(dep_name, '_', ' '), ' Difference between Two Datasets:'],...
                   ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
                   [' - ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_errorbar']));
        end
        
        createFigure;
        plotSimple(indep_vals, difference)  % Plot a simple 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), ' Difference between Two Datasets:'],...
               ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' - ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)

        savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_simple']));
    end

    % Plot 2D data.
    if length(dep_rels1) == 2
        indep_name1 = dep_rels1{1};
        indep_name2 = dep_rels1{2};
        if length(dep_rels2) ~= 2 || ~strcmp(indep_name1, dep_rels2{1}) || ~strcmp(indep_name2, dep_rels2{2})
            error('The selected files do not match.')
        end
        indep_vals1 = data1.(indep_name1);
        indep_vals2 = data1.(indep_name2);
        if ~isfield(data2, indep_name1) ||...
                length(indep_vals1) ~= length(data2.(indep_name1)) ||...
                any(indep_vals1 ~= data2.(indep_name1)) ||...
                ~isfield(data2, indep_name2) ||...
                length(indep_vals2) ~= length(data2.(indep_name2)) ||...
                any(indep_vals2 ~= data2.(indep_name2))
            error('The selected files do not match.')
        end
        
        xunits = getUnits(data1, indep_name1);
        yunits = getUnits(data1, indep_name2);
        zunits = getUnits(data1, dep_name);
        
        difference = dep_vals1 - dep_vals2;
        % Plot the data as a smooth surface.
        if isempty(strfind(indep_name1, 'Phase')) && isempty(strfind(indep_name2, 'Phase'))
            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, difference);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(dep_name, '_', ' '), zunits, ' Difference between Two Datasets:'],...
                   ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
                   [' - ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_smooth']));
        else % Create a polar (smooth) plot.
            if ~isempty(strfind(indep_name1, 'Phase'))
                phase = indep_vals1;
                radius = indep_vals2;
                vals = difference';
                indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                    yunits, '; Phase: ', strrep(indep_name1, '_', ' '), xunits];
            else
                phase = indep_vals2;
                radius = indep_vals1;
                vals = difference;
                indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                    xunits, '; Phase: ', strrep(indep_name2, '_', ' '), yunits];
            end
            createFigure;
            plotPolar(radius, phase, vals);
            title({[strrep(dep_name, '_', ' '), zunits, ' Difference between Two Datasets:'],...
                   ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
                   [' - ', filenames{2}, ' [', data2.Timestamp, ']'], indep_vars}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_smooth']));
        end
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, difference');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), zunits, ' Difference between Two Datasets:'],...
               ['   ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' - ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_pixelated']));
    end
    if length(dep_rels1) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than two sweep variables. ',...
              'The data will not be plotted.'])
    end
end