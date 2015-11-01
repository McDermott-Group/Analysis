function plotQuantumEfficiencyFPT(calibration_coeff, kappa)
% plotQuantumEfficiencyFPT Qunatum efficiency estimation based on the two
% data sets
%   plotQuantumEfficiencyFPT(CALIBRATION_COEFFICIENT, KAPPA) plots quantum
%   efficiency estimated from two text data sets. CALIBRATION_COEFFICIENT
%   [photons/(DAC units)^2] is defined as
%   number_of_photons_inside_the_cavity = CALIBRATION_COEFFICIENT * DAC_amplitude^2
%   at some specific steady-state readout power. KAPPA [sec^-1] is the
%   cavity decay time constant characterizing the energy (photon number) loss.

if ~exist('calibration_coeff', 'var')
    error('Calibration coefficient should be specified as the function argument.')
end

if ~exist('kappa', 'var')
    error('Decay rate kappa should be specified as the function argument.')
end

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

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

if ~isfield(data1, 'Fast_Pulse_Time') || ~isfield(data2, 'Fast_Pulse_Time') ||...
    length(data1.Fast_Pulse_Time) ~= length(data2.Fast_Pulse_Time) ||...
    any(data1.Fast_Pulse_Time ~= data2.Fast_Pulse_Time)
    error('Fast Pulse Time data are not properly specified (missing or unequal).')
end

if ~isfield(data1, 'Readout_Amplitude') || ~isfield(data2, 'Readout_Amplitude') ||...
    length(data1.Readout_Amplitude) ~= length(data2.Readout_Amplitude)
    error('Readout Amplitude data are not properly specified (missing or unequal).')
end

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
        
        xunits = getUnits(data, indep_name);
        
        dep_vals1(dep_vals1 < dep_vals2) = dep_vals2(dep_vals1 < dep_vals2);
        quant_eff = log((1 - dep_vals2) ./ (1 - dep_vals1)) ./...
            (1e-9 * kappa * data1.Fast_Pulse_Time .*...
            (calibration_coeff * data1.Readout_Amplitude.^2));
        quant_eff(~isfinite(quant_eff)) = 0;
        quant_eff(quant_eff < 0) = 0;
        quant_eff(quant_eff > 1) = 1;
        
        if isfield(data1, 'error') && isfield(data1.error, dep_name) &&...
           isfield(data2, 'error') && isfield(data2.error, dep_name) % Plot an errobar graph.
            quant_eff_error = 1.96 * sqrt(data1.error.(dep_name).^2 ./ (1 - dep_vals1).^2 +...
                                          data2.error.(dep_name).^2 ./ (1 - dep_vals2).^2);
                                      
            createFigure('right');
            plotErrorbar(indep_vals, quant_eff, quant_eff_error)
            xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
            ylabel('Quantum Efficiency', 'FontSize', 14);
            title({'Quantum Efficiency Estimated from Two Datasets:',...
                   [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
                   [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
            savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_quant_eff_errorbar']));
        end
        
        createFigure;
        plotSemilog(indep_vals, quant_eff); % Plot a semilog 1D graph.
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel('Quantum Efficiency', 'FontSize', 14);
        title({'Quantum Efficiency Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_quant_eff']));
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
        if strcmp(indep_name1, 'Fast_Pulse_Time')
           quant_eff = log((1 - dep_vals2) ./ (1 - dep_vals1)) ./...
                (1e-9 * kappa * data1.Fast_Pulse_Time(:) * ones(1, size(dep_vals1, 2)) .*...
                (calibration_coeff * data1.Readout_Amplitude.^2));
        elseif strcmp(indep_name2, 'Fast_Pulse_Time')
           quant_eff = log((1 - dep_vals2) ./ (1 - dep_vals1)) ./...
                (1e-9 * kappa * ones(size(dep_vals1, 1), 1) * data1.Fast_Pulse_Time(:)' .*...
                (calibration_coeff * data1.Readout_Amplitude.^2));
        end
        quant_eff(~isfinite(quant_eff)) = 0;
        quant_eff(quant_eff < 0) = 0;
        quant_eff(quant_eff > 1) = 1;
        
        createFigure;
        plotSmooth(indep_vals1, indep_vals2,quant_eff);
        xunits = getUnits(data, indep_name1);
        yunits = getUnits(data, indep_name2);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'Quantum Efficiency Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_quant_eff_smooth']));
        % Plot the data as a pixeleated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, quant_eff');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        title({'Quantum Efficiency Estimated from Two Datasets:',...
               [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
               [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)
        savePlot(fullfile(plts_path, [base_filename1, '_', base_filename2, '_quant_eff_pixelated']));
    end
    if length(dep_rels1) > 2
        disp(['Data variable ''', strrep(dep_name, '_', ' '), ''' depends on more than two sweep variables. ',...
              'The data will not be plotted.'])
    end
end