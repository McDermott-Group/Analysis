function plotFPARFFreqData
%plotFPARFFreqData Show increase in switching probability due to the RF
% excitation. The selected data should be a 2D array of probability values
% with FastPulse Amplitude and RF Frequency as independent variables.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev'))
        continue
    end
    dep_vals = data.(dep_name);
    dep_rels = data.rels.(dep_name);
    
    if isempty(dep_rels)
        disp(['Independent (sweep) variables for data variable ''', strrep(dep_name, '_', ' '), ''' are not specified. ',...
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
            plotErrorbar(indep_vals, dep_vals, data.error.(dep_name));
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title({[filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [filename, '_', dep_name, '_errorbar']));
        end
        
        createFigure;
        plotSimple(indep_vals, dep_vals)  % Plot a simple 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', dep_name, '_delta_prob_simple']));
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        % Plot the data as a smooth surface.
        if ~isempty(strfind(dep_name, 'Probability'))
            dep_vals = dep_vals - mean(dep_vals(:, 1:5), 2) * ones(1, size(dep_vals, 2));
        end
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'RF-Induced Increase in Switching Probability:',...
               [filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', dep_name, '_delta_prob_smooth']));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'RF-Induced Increase in Switching Probability:',...
               [filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', dep_name, '_delta_prob_pixelated']));
    end
    if length(dep_rels) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
end

% Show a message box with the experiment parameters.
showMessageBox(data);