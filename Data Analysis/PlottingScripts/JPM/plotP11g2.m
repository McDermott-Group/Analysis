function plotP11g2
%plotMeasurementData Plot data from a text data file.

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

if isfield(data, 'JPM_A_Switching_Probability') &&...
        isfield(data, 'JPM_B_Switching_Probability')
    PA = data.JPM_A_Switching_Probability;
    PB = data.JPM_B_Switching_Probability;
    PA_PB_prod = PA .* PB;
    if isfield(data.error, 'JPM_A_Switching_Probability') &&...
            isfield(data.error, 'JPM_B_Switching_Probability')
        PA_var = data.error.JPM_A_Switching_Probability.^2;
        PB_var = data.error.JPM_B_Switching_Probability.^2;
        PA_PB_prod_error = sqrt(PA_var .* PB_var +...
            PA_var .* PB.^2 + PA.^2 .* PB_var);
    end
end

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~strcmp(dep_name, 'P11')
        continue
    end
    dep_vals = data.(dep_name);
    dep_rels = data.rels.(dep_name);
    
    if isempty(dep_rels)
        disp(['Independent (sweep) variables for data variable ''', strrep(dep_name, '_', ' '), ''' are not specified.'])
    end

    % Plot 1D P11.
    if length(dep_rels) == 1
        indep_name = dep_rels{1};
        indep_vals = data.(indep_name);
        
        xmin = min(indep_vals);
        xmax = max(indep_vals);
        if xmax == xmin
            xmax = xmax + eps;
        end
        
        xunits = getUnits(data, indep_name);
        yunits = getUnits(data, dep_name);
        
        if isfield(data, 'error') && isfield(data.error, dep_name) % Plot an errobar graph.
            createFigure('right');
            errorbar(indep_vals, dep_vals, 1.96 * data.error.(dep_name),...
                '.', 'LineWidth', 1, 'MarkerSize', 15)
            if exist('PA_PB_prod', 'var')
                hold on
                if exist('PA_PB_prod_error', 'var')
                    errorbar(indep_vals, PA_PB_prod, 1.96 * PA_PB_prod_error,...
                        '.', 'LineWidth', 1, 'MarkerSize', 15)
                else
                    plot(indep_vals, PA_PB_prod,...
                        '.', 'LineWidth', 1, 'MarkerSize', 15)
                end
                hold off
                legend('P_{11}', 'P_A*P_B')
                ymin = min([dep_vals(:) - 1.96 * data.error.(dep_name)(:); PA_PB_prod(:)]);
                ymax = max([dep_vals(:) + 1.96 * data.error.(dep_name)(:); PA_PB_prod(:)]);
            else
                ymin = min([dep_vals(:) - 1.96 * data.error.(dep_name)(:); PA_PB_prod(:)]);
                ymax = max([dep_vals(:) + 1.96 * data.error.(dep_name)(:); PA_PB_prod(:)]);
            end
            if ymin == ymax
                ymax = Inf;
            end
            axis([xmin xmax ymin ymax])
            grid on

            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
            title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_errorbar']));
        end
        
        createFigure;
        plot(indep_vals, dep_vals, '.-', 'LineWidth', 1, 'MarkerSize', 15)  % Plot a simple 1D graph.
        if isfield(data, 'JPM_A_Switching_Probability') && isfield(data, 'JPM_B_Switching_Probability')
            hold on
            plot(indep_vals, PA_PB_prod,...
                '.-', 'LineWidth', 1, 'MarkerSize', 15)
            hold off
            legend('P_{11}', 'P_A*P_B')
            ymin = min([dep_vals(:); PA_PB_prod(:)]);
            ymax = max([dep_vals(:); PA_PB_prod(:)]);
        else
            ymin = min(dep_vals);
            ymax = max(dep_vals);
        end
        if ymin == ymax
            ymax = Inf;
        end
        axis([xmin xmax ymin ymax])
        grid on

        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(dep_name, '_', ' ') yunits], 'FontSize', 14);
        title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_simple']));
    end

    % Plot 2D P11.
    if length(dep_rels) == 2 && strcmp(dep_name, 'P11')
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(indep_name1);
        indep_vals2 = data.(indep_name2);
        
        % Plot the data as a smooth surface.
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_smooth']));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({[strrep(dep_name, '_', ' '), ':'],...
               [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_pixelated']));
    end
    if length(dep_rels) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
              'This data will not be plotted.'])
    end
    
    % Plot 1D g2.
    if isfield(data, 'JPM_A_Switching_Probability') && isfield(data, 'JPM_B_Switching_Probability')
        data.g2 = data.P11 ./ PA_PB_prod;
        
        if length(dep_rels) == 1
            indep_name = dep_rels{1};
            indep_vals = data.(indep_name);
            xunits = getUnits(data, indep_name);

            if isfield(data, 'error') && isfield(data.error, 'g2') % Plot an errobar graph.
                createFigure('right');
                plotErrorbar(indep_vals, data.g2, data.error.g2);
                xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
                ylabel('g_2 = P_{11} / P_A * P_B', 'FontSize', 14);
                title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
                savePlot(fullfile(plts_path, [base_filename, '_g2_errorbar']));
            end

            createFigure;
            plotSimple(indep_vals, data.g2);    % Plot a simple 1D graph.
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel('g_2 = P_{11} / P_A * P_B', 'FontSize', 14);
            title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_g2_simple']));
        end

        % Plot 2D g2.
        if length(dep_rels) == 2 && strcmp(dep_name, 'P11')
            indep_name1 = dep_rels{1};
            indep_name2 = dep_rels{2};
            indep_vals1 = data.(indep_name1);
            indep_vals2 = data.(indep_name2);

            % Plot the data as a smooth surface.
            createFigure;
            plotSmooth(indep_vals1, indep_vals2, data.g2);
            xunits = getUnits(data, indep_name1);
            yunits = getUnits(data, indep_name2);
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({'g2:',...
                   [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_g2_smooth']));
            
            corrected_g2 = data.g2;
            corrected_g2(corrected_g2 > 2) = NaN;
            createFigure('right');
            plotPixelated(indep_vals1, indep_vals2, corrected_g2');
            xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
            ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
            title({'g2:',...
                   [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename, '_g2_pixelated']));
        end
        if length(dep_rels) > 2
            disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
                  'This data will not be plotted.'])
        end
    end
end

% Show a message box with the experiment parameters.
showMessageBox(data);