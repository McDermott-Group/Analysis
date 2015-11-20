function plotFigureOfMerit
%plotFigureOfMerit Plot a figure of merit, specifically, quantum efficiency 
% multiplied by photon flux and divided by the dark tunneling rate. This
% figure could be estimated from two data sets that contain bright
% and dark switching probabilities.

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
    figure_of_merit = log((1 - dep_vals1) ./ (1 - dep_vals2)) ./ log(1 - dep_vals2);
    figure_of_merit(~isfinite(figure_of_merit)) = 0;
    figure_of_merit(figure_of_merit < 0) = 0;
    
    data = data1;
    figure_of_merit_name = 'Figure_of_Merit';
    data.(figure_of_merit_name) = figure_of_merit;
    
    data.units.(figure_of_merit_name) = '';
    data.rels.(figure_of_merit_name) = dep_rels1;
    data.dep{length(data.dep) + 1} = figure_of_merit_name;
    data.plotting.(figure_of_merit_name).full_name =...
        'Figure of Merit (dimensionless): \eta\lambda/\gamma_0';
    data.plotting.(figure_of_merit_name).plot_title =...
        {'Figure of Merit Estimated from Two Datasets:',...
        [' P(bright): ', filenames{1}, ' [', data1.Timestamp, ']'],...
        [' P(dark):   ', filenames{2}, ' [', data2.Timestamp, ']']};
    data.plotting.(figure_of_merit_name).plot_filename =...
            fullfile(plts_path, [base_filename1, '_',...
            base_filename2, '_figofmer']);
    
    plotDataVar(data, figure_of_merit_name)
end