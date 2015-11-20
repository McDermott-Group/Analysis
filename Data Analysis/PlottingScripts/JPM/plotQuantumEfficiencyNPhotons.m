function plotQuantumEfficiencyNPhotons(number_of_photons)
%plotQuantumEfficiencyNPhotons(NUMBER_OF_PHOTONS) Plot the quantum efficiency 
% estimated from two text data sets based on specified number of photons
% in the measuremnet interval (NUMBER_OF_PHOTONS).

if ~exist('number_of_photons', 'var')
    error('Number of photons in the measurment interval should be specified as the function argument.')
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
     
    if isempty(dep_rels1) || isempty (dep_rels2)
        continue
    end
    
    if length(dep_rels1) ~= length(dep_rels2)
        error('The selected files do not match.')
    end
    for k = 1:length(dep_rels1)
        if ~strcmp(dep_rels1{k}, dep_rels2{k}) ||...
               length(dep_rels1{k}) ~= length(dep_rels2{k}) ||...
               any(dep_rels1{k} ~= dep_rels2{k})
            error('The selected files do not match.')
        end
    end

    dep_vals1(dep_vals1 < dep_vals2) = dep_vals2(dep_vals1 < dep_vals2);
    quant_eff = log((1 - dep_vals2) ./ (1 - dep_vals1)) / number_of_photons;
    quant_eff(dep_vals1 == 1) = 0;
    quant_eff(quant_eff < 0) = 0;
    quant_eff(quant_eff > 1) = 1;
    
    data = data1;
    quant_eff_name = 'Quantum_Efficiency';
    data.(quant_eff_name) = quant_eff;
    
    if isfield(data1, 'error') && isfield(data1.error, dep_name) &&...
       isfield(data2, 'error') && isfield(data2.error, dep_name)
        data.error.(quant_eff_name) =...
                sqrt(data1.error.(dep_name).^2 ./ (1 - dep_vals1).^2 +...
                     data2.error.(dep_name).^2 ./ (1 - dep_vals2).^2);
    end
    
    data.units.(quant_eff_name) = '';
    data.rels.(quant_eff_name) = dep_rels1;
    data.dep{length(data.dep) + 1} = quant_eff_name;
    data.plotting.(quant_eff_name).plot_title =...
        {'Quantum Efficiency Estimated from Two Datasets:',...
        [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
        [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']};
    data.plotting.(quant_eff_name).plot_filename =...
            fullfile(plts_path, [base_filename1, '_',...
            base_filename2, '_', quant_eff_name]);
    
    plotDataVar(data, quant_eff_name)
end