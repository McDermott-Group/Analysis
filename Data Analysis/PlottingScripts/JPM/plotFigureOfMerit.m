function plotFigureOfMerit
%plotFigureOfMerit Plot a figure of merit, specifically, quantum efficiency 
% multiplied by photon flux and divided by the dark tunneling rate. This
% figure could be estimated from two text data sets: data sets with bright
% and dark switching probabilities.

% Select files for quantum efficiency estimation.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file containing P[bright]...',...
     'Select the second data file containing P[dark]...'});
if ~status
    return
end

% Read probability data file, convert the variable names, and define
% the units.
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
    if ~strcmp(dep_name, 'Switching_Probability') &&...
            ~strcmp(dep_name, 'JPM_A_Switching_Probability') &&...
            ~strcmp(dep_name, 'JPM_B_Switching_Probability')
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
        disp(['Independent (sweep) variables for data variable ''', strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end
    
    % Plot 1D data.
    if length(dep_rels1) == 1
        indep_name = dep_rels1{1};
        if length(dep_rels2) ~= 1 || ~strcmp(indep_name, dep_rels2{1})
            error('The selected files do not match.')
        end
        indep_vals = data1.(indep_name);
        if ~isfield(data2, indep_name)  || length(indep_vals) ~= length(data2.(indep_name)) ||...
                any(indep_vals ~= data2.(indep_name))
            error('The selected files do not match.')
        end

        xunits = getUnits(data1, indep_name);
        
        figure_of_merit = log((1 - dep_vals1) ./ (1 - dep_vals2)) ./ log(1 - dep_vals2);
        figure_of_merit(dep_vals1 == 1) = 0;
        figure_of_merit(figure_of_merit < 0) = 0;
        
        createFigure;
        plotSemilog(indep_vals, figure_of_merit);  % Plot a semilog 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel('Figure of Merit: \eta\lambda/\gamma_0  (dimensionless)', 'FontSize', 14);
        title({'Figure of Merit Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_fom']));
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
        
        % Plot the data as a smooth surface.
        dep_vals1(dep_vals1 < dep_vals2) = dep_vals2(dep_vals1 < dep_vals2);
        figure_of_merit = log((1 - dep_vals1) ./ (1 - dep_vals2)) ./ log(1 - dep_vals2);
        figure_of_merit(dep_vals1 == 1) = 0;
        figure_of_merit(figure_of_merit < 0) = 0;
     
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, figure_of_merit);
     
        xunits = getUnits(data1, indep_name1);
        yunits = getUnits(data1, indep_name2);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'Figure of Merit Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_fom_smooth']));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, figure_of_merit');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'Figure of Merit Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_fom_pixelated']));
    end
    if length(dep_rels1) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
              'The data will not be plotted.'])
    end
end