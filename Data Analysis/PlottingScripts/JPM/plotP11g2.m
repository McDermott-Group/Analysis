function plotP11g2
%plotMeasurementData Plot data from a text data file.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

if ~isfield(data, 'P11')
    error('The selected data set does not contain joint switching probability.')
end
if ~isfield(data, 'JPM_A_Switching_Probability') &&...
        isfield(data, 'JPM_B_Switching_Probability')
    error('The selected data set does not contain individual switching probabilities.')
end

P11 = data.P11;
dep_rels = data.rels.P11;
if isfield(data.error, 'P11')
    P11_var = data.error.P11.^2;
end


PA = data.JPM_A_Switching_Probability;
PB = data.JPM_B_Switching_Probability;
PA_PB_prod = PA .* PB;
g2 = P11 ./ PA_PB_prod;

if isfield(data.error, 'JPM_A_Switching_Probability') &&...
        isfield(data.error, 'JPM_B_Switching_Probability')
    PA_var = data.error.JPM_A_Switching_Probability.^2;
    PB_var = data.error.JPM_B_Switching_Probability.^2;
    PA_PB_prod_var = (PA_var .* PB_var +...
        PA_var .* PB.^2 + PA.^2 .* PB_var);
    if exist('P11_var', 'var')
        g2_var = g2.^2 .* (P11_var ./ P11.^2 + PA_PB_prod_var ./ PA_PB_prod.^2);
    end
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
    yunits = getUnits(data, 'P11');

    if isfield(data, 'error') && isfield(data.error, 'P11') % Plot an errobar graph.
        createFigure('right');
        errorbar(indep_vals, P11, 1.96 * data.error.P11,...
            '.', 'LineWidth', 1, 'MarkerSize', 15)
        if exist('PA_PB_prod', 'var')
            hold on
            if exist('PA_PB_prod_var', 'var')
                PA_PB_prod_error = 1.96 * sqrt(PA_PB_prod_var);
                errorbar(indep_vals, PA_PB_prod, PA_PB_prod_error,...
                    '.', 'LineWidth', 1, 'MarkerSize', 15)
            else
                plot(indep_vals, PA_PB_prod,...
                    '.', 'LineWidth', 1, 'MarkerSize', 15)
                PA_PB_prod_error = zeros(size(PA_PB_prod));
            end
            hold off
            legend('P_{11}', 'P_A*P_B')
            ymin = min([P11(:) - 1.96 * data.error.P11(:);...
                PA_PB_prod(:) - PA_PB_prod_error(:)]);
            ymax = max([P11(:) + 1.96 * data.error.P11(:);...
                PA_PB_prod(:) - PA_PB_prod_error(:)]);
        else
            ymin = min(P11(:) - 1.96 * data.error.P11(:));
            ymax = max(P11(:) + 1.96 * data.error.P11(:));
        end
        if ymin == ymax
            ymax = Inf;
        end
        axis([xmin xmax ymin ymax])
        grid on

        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel('P11', 'FontSize', 14);
        title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', 'P11', '_errorbar']));
    end

    createFigure;
    plot(indep_vals, P11, '.-', 'LineWidth', 1, 'MarkerSize', 15)  % Plot a simple 1D graph.
    if isfield(data, 'JPM_A_Switching_Probability') && isfield(data, 'JPM_B_Switching_Probability')
        hold on
        plot(indep_vals, PA_PB_prod,...
            '.-', 'LineWidth', 1, 'MarkerSize', 15)
        hold off
        legend('P_{11}', 'P_A*P_B')
        ymin = min([P11(:); PA_PB_prod(:)]);
        ymax = max([P11(:); PA_PB_prod(:)]);
    else
        ymin = min(P11);
        ymax = max(P11);
    end
    if ymin == ymax
        ymax = Inf;
    end
    axis([xmin xmax ymin ymax])
    grid on

    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel(['P_{11}', yunits], 'FontSize', 14);
    title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_P11_simple']));
end

% Plot 2D P11.
if length(dep_rels) == 2
    indep_name1 = dep_rels{1};
    indep_name2 = dep_rels{2};
    indep_vals1 = data.(indep_name1);
    indep_vals2 = data.(indep_name2);

    % Plot the data as a smooth surface.
    createFigure;
    plotSmooth(indep_vals1, indep_vals2, P11);
    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
    title({'P11:', [filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_P11_smooth']));
    % Plot the data as a pixeleated image.
    createFigure('right');
    plotPixelated(indep_vals1, indep_vals2, P11');
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
    title({'P11:',...
           [filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_P11_pixelated']));
end
if length(dep_rels) > 2
    disp(['Data variable ''''P11'''' depends on more than two sweep variables. ',...
          'This data will not be plotted.'])
end

% Plot 1D g2.
if length(dep_rels) == 1
    indep_name = dep_rels{1};
    indep_vals = data.(indep_name);
    xunits = getUnits(data, indep_name);

    if exist('g2_var', 'var')
        createFigure('right');
        plotErrorbar(indep_vals, g2, 1.96 * sqrt(g2_var));
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel('g_2 = P_{11} / P_AP_B', 'FontSize', 14);
        title({[filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_g2_errorbar']));
    end

    createFigure;
    plotSimple(indep_vals, g2);    % Plot a simple 1D graph.
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel('g_2 = P_{11} / P_AP_B', 'FontSize', 14);
    title({[filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_g2_simple']));
end

% Plot 2D g2.
if length(dep_rels) == 2
    indep_name1 = dep_rels{1};
    indep_name2 = dep_rels{2};
    indep_vals1 = data.(indep_name1);
    indep_vals2 = data.(indep_name2);

    % Plot the data as a smooth surface.
    createFigure;
    plotSmooth(indep_vals1, indep_vals2, g2);
    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
    title({'g2:',...
           [filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_g2_smooth']));

    corrected_g2 = g2;
    corrected_g2(corrected_g2 > 2) = NaN;
    createFigure('right');
    plotPixelated(indep_vals1, indep_vals2, corrected_g2');
    xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
    title({'g2:',...
           [filename, ext, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_g2_pixelated']));
end

if length(dep_rels) > 2
    disp(['Data variable ''''P11'''' depends on more than two sweep variables.',...
          'Plotting is not implemented.'])
end

% Show a message box with the experiment parameters.
showMessageBox(data);

disp(['<g_2> = ', num2str(mean(g2(:))), ' ± ',...
    num2str(1.96 * std(g2(:)) / sqrt(length(g2(:)))),...
    ' (mean ± 1.96 * standard error of the mean)'])