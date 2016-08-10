function plotMaxIQContrast
%plotMaxIQContrast  Plot the data in the frame that maximizes
% the IQ-space variations.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Maximize IQ-contrast.
data = maximizeIQContrast(data);

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    
    if ~isempty(strfind(I_name, '_Std_Dev')) ||...
            isempty(strfind(I_name, 'I'))
        continue
    else
        Q_name = strrep(I_name, 'I', 'Q');
        dep_name = strrep(I_name, 'I', 'Phase_Space_Shift');
        res_name = strrep(I_name, 'I', 'Residual');
        if ~isfield(data, Q_name) || ~isfield(data, dep_name) ||...
                ~isfield(data, res_name)
            continue
        end
    end

    dep_rels = data.rels.(dep_name);

    dep_vals = data.(dep_name);
    res_vals = data.(res_name);

    % Plot 1D data.
    if length(dep_rels) == 1
        indep_name = dep_rels{1};
        indep_vals = data.(indep_name);
        
        xunits = getUnits(data, indep_name);
        yunits = getUnits(data, I_name);
        
        % Plot an errobar graph.
        if isfield(data, 'error') && isfield(data.error, I_name) &&...
                isfield(data.error, Q_name)
            dep_err = 1.96 * data.error.(dep_name);
            res_err = 1.96 * data.error.(res_name);
            createFigure('right');
            plotErrorbar(indep_vals, res_vals, res_err)
            hold on
                 plotErrorbar(indep_vals, dep_vals, dep_err)
            hold off
            legend(strrep(res_name, '_', ' '), strrep(dep_name, '_', ' '))
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
            ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
            title({[strrep(filename, '_', '\_'), ext,...
                ' [', data.Timestamp, ']']}, 'FontSize', 10)
            savePlot(fullfile(plts_path, [filename, '_', dep_name,...
                '_rotframe_errorbar']));
        end
        
        createFigure;
        plotSimple(indep_vals, res_vals)  % Plot a simple 1D graph.
        hold on
            plotSimple(indep_vals, dep_vals)
        hold off
        legend(strrep(res_name, '_', ' '), strrep(dep_name, '_', ' '))
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(dep_name, '_', ' '), yunits], 'FontSize', 14)
        title({[filename, ext, ' [', data.Timestamp, ']']},...
                'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', dep_name,...
            '_rotframe_simple']));
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        zunits = getUnits(data, I_name);

        if isempty(strfind(indep_name1, 'Phase')) &&...
                isempty(strfind(indep_name2, 'Phase'))
            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, dep_vals);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
            title({[strrep(dep_name, '_', ' '), zunits, ':'],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
            savePlot(fullfile(plts_path, [filename, '_', dep_name,...
                '_rotframe_smooth']));
            
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, res_vals);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(res_name, '_', ' '), zunits, ':'],...
                   [filename, ext, ' [', data.Timestamp, ']']},...
                   'Interpreter', 'none', 'FontSize', 10)
        else % Create a polar (smooth) plot.
            if ~isempty(strfind(indep_name1, 'Phase'))
                phase = indep_vals1;
                radius = indep_vals2;
                vals = dep_vals';
                res_vals = res_vals';
                indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                    yunits, '; Phase: ', strrep(indep_name1, '_', ' '),...
                    xunits];
            else
                phase = indep_vals2;
                radius = indep_vals1;
                vals = dep_vals;
                indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                    xunits, '; Phase: ', strrep(indep_name2, '_', ' '),...
                    yunits];
            end
            createFigure;
            plotPolar(radius, phase, vals);
            title({[strrep(dep_name, '_', ' '), zunits, ':'],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']'], indep_vars}, 'FontSize', 10)
            savePlot(fullfile(plts_path, [filename, '_', dep_name,...
                '_rotframe_smooth']));
            
            createFigure;
            plotPolar(radius, phase, res_vals);
            title({[res_name, zunits, ':'],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']'], indep_vars},...
                   'FontSize', 10)
        end
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        title({[strrep(dep_name, '_', ' '), zunits, ':'],...
               [strrep(filename, '_', '\_'), ext,...
               ' [', data.Timestamp, ']']}, 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', dep_name,...
            '_rotframe_pixelated']));
        
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, res_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        title({[strrep(res_name, '_', ' '), zunits, ':'],...
               [strrep(filename, '_', '\_'), ext,...
               ' [', data.Timestamp, ']']}, 'FontSize', 10)
    end
end
end