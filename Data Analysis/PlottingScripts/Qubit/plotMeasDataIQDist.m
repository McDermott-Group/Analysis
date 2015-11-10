function plotMeasDataIQDist
%plotMeasDataIQDist Plot IQ-distance between two data sets.

% Select files to compute the difference.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file...',...
     'Select the second data file...'});
if ~status
    return
end

% Read the data files, convert the variable names, and specify the units.
data1 = processMeasurementData(importMeasurementData(fullfile(pathnames{1}, filenames{1})));
data2 = processMeasurementData(importMeasurementData(fullfile(pathnames{2}, filenames{2})));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

for data_index = 1:length(data1.dep)
    I_name = data1.dep{data_index};
    if ~isempty(strfind(I_name, '_Std_Dev'))
        continue
    end
    
    if ~isempty(strfind(I_name, 'I'))
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data1, Q_name)
            continue
        end
        dep_name = strrep(I_name, 'I', 'IQ-Distance');
    else
        continue
    end

    if ~isfield(data2, I_name) && ~isfield(data2, Q_name)
        error('The selected files do not match.')
    else
        I1 = data1.(I_name);
        Q1 = data1.(Q_name);
        I2 = data2.(I_name);
        Q2 = data2.(Q_name);
    end
    I1_rels = data1.rels.(I_name);
    Q1_rels = data1.rels.(Q_name);
    I2_rels = data2.rels.(I_name);
    Q2_rels = data2.rels.(Q_name);
    
    if isempty(I1_rels) || isempty(Q1_rels) ||...
        isempty(I2_rels) || isempty(Q2_rels)
        continue
    end
    
    dim = length(I1_rels);
    if length(I2_rels) ~= dim ||...
            length(Q1_rels) ~= dim || length(Q2_rels) ~= dim
        error('The selected files do not match.')
    end
    for k = 1:length(I1_rels)
        if ~strcmp(I1_rels{k}, I2_rels{k}) ||...
               ~strcmp(Q1_rels{k}, Q2_rels{k}) ||...
               ~strcmp(I1_rels{k}, Q1_rels{k})
            error('The selected files do not match.')
        end
    end
    
    try
        distance = sqrt((I2 - I1).^2 + (Q2 - Q1).^2);
    catch
        error('The selected files do not match.')
    end
    
    % Plot 1D data.
    if dim == 1
        indep_name = I1_rels{1};
        indep_vals = data1.(indep_name);

        xunits = getUnits(data1, indep_name);
        yunits = getUnits(data1, I_name);
        
        createFigure;
        plotSimple(indep_vals, distance)  % Plot a simple 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), ' between Two Datasets:'],...
               [filenames{1}, ' [', data1.Timestamp, ']'],...
               [filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)

        savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name]));
    end

    % Plot 2D data.
    if dim == 2
        indep_name1 = I1_rels{1};
        indep_name2 = I1_rels{2};

        indep_vals1 = data1.(indep_name1);
        indep_vals2 = data1.(indep_name2);
        
        xunits = getUnits(data1, indep_name1);
        yunits = getUnits(data1, indep_name2);
        zunits = getUnits(data1, dep_name);
        
        % Plot the data as a smooth surface.
        if isempty(strfind(indep_name1, 'Phase')) && isempty(strfind(indep_name2, 'Phase'))
            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, distance);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(dep_name, '_', ' '), zunits, ' between Two Datasets:'],...
                   [filenames{1}, ' [', data1.Timestamp, ']'],...
                   [filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_smooth']));
        else % Create a polar (smooth) plot.
            if ~isempty(strfind(indep_name1, 'Phase'))
                phase = indep_vals1;
                radius = indep_vals2;
                vals = distance';
                indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                    yunits, '; Phase: ', strrep(indep_name1, '_', ' '), xunits];
            else
                phase = indep_vals2;
                radius = indep_vals1;
                vals = distance;
                indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                    xunits, '; Phase: ', strrep(indep_name2, '_', ' '), yunits];
            end
            createFigure;
            plotPolar(radius, phase, vals);
            title({[strrep(dep_name, '_', ' '), zunits, ' between Two Datasets:'],...
                   [filenames{1}, ' [', data1.Timestamp, ']'],...
                   [filenames{2}, ' [', data2.Timestamp, ']'], indep_vars}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_polar']));
        end
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, distance');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), zunits, ' between Two Datasets:'],...
               [filenames{1}, ' [', data1.Timestamp, ']'],...
               [filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '-', base_filename2, '_', dep_name, '_pixelated']));
    end
    if length(dim) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than two sweep variables. ',...
              'The data will not be plotted.'])
    end
end