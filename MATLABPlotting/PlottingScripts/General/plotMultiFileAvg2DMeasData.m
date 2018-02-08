function plotMultiFileAvg2DMeasData(data_variable)
% plotMultiFileAvg2DMeasData(DATA_VARIABLE) Plots a pixel-by-pixel mean of
%   several 2D sweep data files
% 
% 



% Select files to plot.
[filenames, pathnames, status] = selectMeasurementDataFile;
if ~status
    return
end

if ~iscell(filenames)
    filenames = {filenames};
end
if ~iscell(pathnames)
    pathnames = {pathnames};
end

% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end

if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
elseif ~iscell(data_variable)
    dep_vars = {data_variable};
else
    dep_vars = data_variable;
end
shapes=size(data{1}.(dep_vars{1}));
dep_vals = zeros(length(filenames),shapes(1),shapes(2));
for k = 1:length(dep_vars)
    data_variable = dep_vars{k};
    
    % assume that the dep_vars for each file are the same
    dep_rels1 = data{1}.rels.(data_variable);
    
    for q = 1:length(filenames)
        % Check that the data variable exists (compute it if necessary).
        [data{q}, data_variable] = checkDataVar(data{q}, data_variable);
        dep_vals(q,:,:) = data{q}.(data_variable);
        % dep_rels1 = data{q}.rels.(data_variable);
        dep_rels = data{q}.rels.(data_variable);
        
        if ~isequal(dep_rels,dep_rels1)
            error('Independent variables are not the same for all selected files.');
        end
        
        % data{1}.rels.(data_variable)(:,1) = [];


%         if ~isfield(data{q}, data_variable)
%             error(['Could not find data variable ''', data_variable, ''' in ',...
%                 'at least one of the files.'])
%         elseif isempty(dep_rels)
%                 error(['Independent (sweep) variables for data variable ''',...
%                       strrep(data_variable, '_', ' '),...
%                       ''' are not specified.'])
%         end
    end
    % dep_vals = mean(dep_vals,1);
    avg_data = squeeze(mean(dep_vals,1));

    
    % avg_data(:,:) = avg_data1(:,:,:);
    processed_data_var = ['Average_', data_variable];
    data2plot = data{1};
    data2plot.(processed_data_var) = avg_data;
    data2plot.units.(processed_data_var) = data{1}.units.(data_variable);
    data2plot.rels.(processed_data_var) = data{1}.rels.(data_variable);
    data2plot.dep{length(data{1}.dep) + 1} = processed_data_var;
    units = getUnits(data{1}, data_variable);
    data2plot.plotting.(processed_data_var).plot_title =...
            {[strrep(dep_vars{k}, '_', ' '), units],...
             [strrep(filenames{1}, '_', '\_'), ' - ',...
              strrep(filenames{end}, '_', '\_')],...
             ['[', data{1}.Timestamp, ' - ',...
              data{end}.Timestamp, ']']};
    data2plot.plotting.(processed_data_var).full_name =...
        ['_'];
    data2plot.plotting.(processed_data_var).extra_filename =...
        ['_'];

    plotDataVar(data2plot, processed_data_var);
    
    saveMeasData(data2plot, [data_variable, '_avg'])
end
end