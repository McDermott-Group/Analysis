function plotMeasData
%plotMeasData Plot data from a text data file.

% Select a file to plot.
[filename, pathname, status] = selectMeasurementDataFile(1);
if ~status
    return
end

% Read the data file, convert the variable names, and specify the units.
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
    
    if isempty(dep_rels)
        disp(['Independent (sweep) variables for data variable ''',...
              strrep(dep_name, '_', ' '), ''' are not specified. ',...
              'This data will not be plotted.'])
    end
    
    % Plot 1D data.
    if length(dep_rels) == 1
        indep_name = dep_rels{1};
        indep_vals = data.(indep_name);
        
        xunits = getUnits(data, indep_name);
        yunits = getUnits(data, dep_name);
        
        if isfield(data, 'error') && isfield(data.error, dep_name) % Plot an errobar graph.
            createFigure('right');
            plotErrorbar(indep_vals, dep_vals, data.error.(dep_name))
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_errorbar']));
        end
        
        createFigure;
        plotSimple(indep_vals, dep_vals)  % Plot a simple 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_simple']));
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        zunits = getUnits(data, dep_name);

        if isempty(strfind(indep_name1, 'Phase')) && isempty(strfind(indep_name2, 'Phase'))
            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, dep_vals);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(dep_name, '_', ' '), zunits, ':'],...
                   [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_smooth']));
        else % Create a polar (smooth) plot.
            if ~isempty(strfind(indep_name1, 'Phase'))
                phase = indep_vals1;
                radius = indep_vals2;
                vals = dep_vals';
                indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                    yunits, '; Phase: ', strrep(indep_name1, '_', ' '), xunits];
            else
                phase = indep_vals2;
                radius = indep_vals1;
                vals = dep_vals;
                indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                    xunits, '; Phase: ', strrep(indep_name2, '_', ' '), yunits];
            end
            createFigure;
            plotPolar(radius, phase, vals);
            title({[strrep(dep_name, '_', ' '), zunits, ':'],...
                   [filename, ' [', data.Timestamp, ']'], indep_vars},...
                   'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_smooth']));
        end
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), zunits, ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_pixelated']));
    end
    if length(dep_rels) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
end

% Show a message box with the experiment parameters.
showMessageBox(data);