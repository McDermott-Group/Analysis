function plotDoubleJunctionJPMCritCurrent(I_critical, R)
%plotDoubleJunctionJPMCritCurrent Plot critical current vs differential
%bias
%   plotDoubleJunctionJPMCritCurrent(I_CRITICAL, R) plots critical current
%   vs differential bias. I_critical is the critical current in Amps of
%   a single junction. If no critical current is passed to the function
%   the data will be plotted against measured values. R is the current 
%   bias resistance in Ohms (1000 Ohm by default). 

plot_meas = false;

if ~exist('I_critical', 'var')
    plot_meas = true;
end

if ~exist('R', 'var')
    R = 1e3; % Ohm 
end

% Select files for the plot.
[filename, pathname, status] = selectMeasurementDataFile(1,...
    {'Select a data file with the Input and Output Bias Voltages as independent axes...'});
if ~status
    return
end

% Read the data file, convert the variable names, and specify the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

[~, base_filename] = fileparts(filename);

if ~isfield(data, 'Bias_Voltage') && ~isempty(data.Bias_Voltage)
    error('Bias Voltage data is missing.')
end

if ~isfield(data, 'Input_Bias_Voltage') && ~isempty(data.Input_Bias_Voltage)
    error('Input Bias Voltage data is missing.')
end

if isfield(data, 'Switching Probability')
    data_name = 'Switching Probability';
elseif isfield(data, 'P11')
    data_name = 'P11';
else
    error(['Could not find "Switching Probability" or "P11" variables' ...
        'in data file.'])
end    

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~strcmp(dep_name, data_name)
        continue
    end
    dep_rels = data.rels.(dep_name);

    [rows, cols] = find(data.(dep_name) == 1);
    if strcmp(dep_rels{1}, 'Bias_Voltage')
        input_bias = data.Input_Bias_Voltage(cols) / R;
        output_bias = data.Bias_Voltage(rows) / R;
    elseif strcmp(dep_rels{1}, 'Input_Bias_Voltage')
        input_bias = data.Input_Bias_Voltage(rows) / R;
        output_bias = data.Bias_Voltage(cols) / R;
    end
    
    createFigure;
    if plot_meas
        plotDots(output_bias*R, input_bias*R);
        xlabel('Voltage Bias (V)', 'FontSize', 14);
        ylabel('Input Voltage Bias (V)', 'FontSize', 14);
    else 
        plotDots((output_bias - input_bias) / (2 * I_critical),...
            (input_bias + output_bias) / I_critical)
        ylim([0 2.5])
        xlabel('\Delta=(I_{output}-I_{input})/(2I_{critical})', ...
            'FontSize', 14);
        ylabel('(I_{output}+I_{input})/I_{critical}', 'FontSize', 14);
    end
    title([filename, ' [', data.Timestamp, ']'], 'Interpreter', 'none', 'FontSize', 10)
    savePlot(fullfile(plts_path, [base_filename, '_', dep_name, '_pixelated']));
end