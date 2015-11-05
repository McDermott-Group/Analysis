function plotMaxIQContrast
%plotMaxIQContrast Plot the data in the frame that maximizes variations.

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
    I_name = data.dep{data_index};
    if ~isempty(strfind(I_name, '_Std_Dev'))
        continue
    end
    
    if ~isempty(strfind(I_name, 'I'))
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data, Q_name)
            continue
        end
        dep_name = strrep(I_name, 'I', 'Maximum IQ-Contrast');
        res_name = 'Residual';
    else
        continue
    end

    I = data.(I_name);
    Q = data.(Q_name);

    I_rels = data.rels.(I_name);
    Q_rels = data.rels.(Q_name);
    
    if isempty(I_rels) || isempty(Q_rels)
        continue
    end
    
    I = I - mean(I(:));
    Q = Q - mean(Q(:));
    phi = findAngle(I(:), Q(:));
    dep_vals = I * cos(phi) - Q * sin(phi);
    res_vals = I * sin(phi) + Q * cos(phi);
    
    % Plot 1D data.
    if length(I_rels) == 1 && length(Q_rels) == 1
        indep_name = I_rels{1};
        indep_vals = data.(indep_name);
        
        xunits = getUnits(data, indep_name);
        yunits = getUnits(data, I_name);
        
        % Plot an errobar graph.
        if isfield(data, 'error') && isfield(data.error, I_name) &&...
                isfield(data.error, Q_name)
            dep_err = 1.96 * sqrt(data.error.(I_name).^2 * cos(phi)^2 +...
                                data.error.(Q_name).^2 * sin(phi)^2);
            res_err = 1.96 * sqrt(data.error.(I_name).^2 * sin(phi)^2 +...
                    data.error.(Q_name).^2 * cos(phi)^2);
            createFigure('right');
            plotErrorbar(indep_vals, res_vals, res_err)
            hold on
                 plotErrorbar(indep_vals, dep_vals, res_err)
            hold off
            legend(res_name, dep_name)
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_rotframe_errorbar']));
        end
        
        createFigure;
        plotSimple(indep_vals, res_vals)  % Plot a simple 1D graph.
        hold on
            plotSimple(indep_vals,  dep_vals)
        hold off
        legend(res_name, dep_name)
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_rotframe_simple']));
    end

    % Plot 2D data.
    if length(I_rels) == 2 && length(Q_rels) == 2
        indep_name1 = I_rels{1};
        indep_name2 = Q_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        zunits = getUnits(data, I_name);

        if isempty(strfind(indep_name1, 'Phase')) && isempty(strfind(indep_name2, 'Phase'))
            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, dep_vals);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(dep_name, '_', ' '), zunits, ':'],...
                   [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_rotframe_smooth']));
            
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, res_vals);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({[strrep(res_name, '_', ' '), zunits, ':'],...
                   [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        else % Create a polar (smooth) plot.
            if ~isempty(strfind(indep_name1, 'Phase'))
                phase = indep_vals1;
                radius = indep_vals2;
                vals = dep_vals';
                res_vals = res_vals';
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
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_rotframe_smooth']));
            
            createFigure;
            plotPolar(radius, phase, res_vals);
            title({[res_name, zunits, ':'],...
                   [filename, ' [', data.Timestamp, ']'], indep_vars},...
                   'Interpreter', 'none', 'FontSize', 10)
        end
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), zunits, ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_rotframe_pixelated']));
        
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, res_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(res_name, '_', ' '), zunits, ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    end
    if length(I_rels) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '),...
              ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
end
end

function phi = findAngle(I, Q)

phi = fminsearch(@costFun, 0);

    function SumSquaredError  = costFun(phi)
        SumSquaredError = sum((I * sin(phi) + Q * cos(phi)).^2);
    end
end