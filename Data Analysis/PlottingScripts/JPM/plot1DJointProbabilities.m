function plot1DJointProbabilities
%plot1DJointProbabilies Plot probabilities from a data file.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Sanity checks.
if ~isfield(data, 'JPM_A_Switching_Probability') ||...
        length(data.rels.JPM_A_Switching_Probability) ~= 1 ||...
        ~isfield(data, 'JPM_B_Switching_Probability') ||...
        length(data.rels.JPM_B_Switching_Probability) ~= 1 ||...
        ~isfield(data, 'P00') || length(data.rels.P00) ~= 1 ||...
        ~isfield(data, 'P01') || length(data.rels.P01) ~= 1 ||...
        ~isfield(data, 'P10') || length(data.rels.P10) ~= 1 ||...
        ~isfield(data, 'P11') || length(data.rels.P11) ~= 1
    error(['The data should have the following 1D data variables ',...
        '''JPM_A_Switching_Probability'', ''JPM_B_Switching_Probability''',...
        ''', ''P00'', ''P01'', ''P10'', and ''P11''. These data ',...
        'variables should depend on an single independent (sweep) variable.'])
end

if ~strcmp(data.rels.P00{1}, data.rels.JPM_A_Switching_Probability{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.JPM_B_Switching_Probability{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P01{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P10{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P11{1})
    error(['The independent (sweep) variable should be the same for all',...
        ' probability data variables.'])
end

createFigure;

plot(data.(data.rels.JPM_A_Switching_Probability{1}),...
     data.JPM_A_Switching_Probability, 'r.-',...
     data.(data.rels.JPM_B_Switching_Probability{1}),...
     data.JPM_B_Switching_Probability, 'b.-',...
     data.(data.rels.P00{1}), data.P00, 'k.-',...
     data.(data.rels.P01{1}), data.P01, 'g.-',...
     data.(data.rels.P10{1}), data.P10, 'c.-',...
     data.(data.rels.P11{1}), data.P11, 'm.-',...
     'LineWidth', 1, 'MarkerSize', 15)

ylim([0 1])
grid on
legend('P_{JPM A}', 'P_{JPM B}',...
       'P_{00}', 'P_{01}', 'P_{10}',...
       'P_{11}')
xunits = getUnits(data, data.rels.P00{1});
xlabel([strrep(data.rels.P00{1}, '_', ' '), xunits], 'FontSize', 14)
ylabel('Switching Probability', 'FontSize', 14)
title({[filename, ext, ' [', data.Timestamp, ']']},...
    'Interpreter', 'none', 'FontSize', 10)

savePlot(fullfile(plts_path, [filename, '_prob_simple']));

createFigure('right');

hold on
if isfield(data.error, 'JPM_A_Switching_Probability')
    errorbar(data.(data.rels.JPM_A_Switching_Probability{1}),...
        data.JPM_A_Switching_Probability,...
        1.96 * data.error.JPM_A_Switching_Probability,...
        'r.', 'LineWidth', 1, 'MarkerSize', 15)
end
if isfield(data.error, 'JPM_B_Switching_Probability')
    errorbar(data.(data.rels.JPM_B_Switching_Probability{1}),...
        data.JPM_B_Switching_Probability,...
        1.96 * data.error.JPM_B_Switching_Probability,...
        'b.', 'LineWidth', 1, 'MarkerSize', 15)
end
if isfield(data.error, 'P00')
    errorbar(data.(data.rels.P00{1}), data.P00,...
        1.96 * data.error.P00, 'k.', 'LineWidth', 1, 'MarkerSize', 15)
end
if isfield(data.error, 'P01')
    errorbar(data.(data.rels.P01{1}), data.P01,...
        1.96 * data.error.P01, 'g.', 'LineWidth', 1, 'MarkerSize', 15)
end
if isfield(data.error, 'P10')
    errorbar(data.(data.rels.P10{1}), data.P10,...
        1.96 * data.error.P10, 'c.', 'LineWidth', 1, 'MarkerSize', 15)
end
if isfield(data.error, 'P11')
    errorbar(data.(data.rels.P11{1}), data.P11,...
        1.96 * data.error.P11, 'm.', 'LineWidth', 1, 'MarkerSize', 15)
end
hold off

ylim([0 1])
grid on
legend('P_{JPM A}', 'P_{JPM B}',...
       'P_{00}', 'P_{01}', 'P_{10}',...
       'P_{11}')
xunits = getUnits(data, data.rels.P00{1});
xlabel([strrep(data.rels.P00{1}, '_', ' '), xunits], 'FontSize', 14)
ylabel('Switching Probability', 'FontSize', 14)
title({[filename, ext, ' [', data.Timestamp, ']']},...
    'Interpreter', 'none', 'FontSize', 10)

savePlot(fullfile(plts_path, [filename, '_prob_errobar']));

% Show a message box with the experiment parameters.
showMessageBox(data);