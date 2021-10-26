function plotMultiFileAvg1DMeasData(data_variable)
% plotMultiFileAvg2DMeasData(
% DATA_VARIABLE) Plots a mean of several 1D data files 



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

% indep_name = data{1}.rels.(dep_name){1}
% for k = 1:length(filenames)
%     if  ~strcmp(indep_name, data{k}.rels.(dep_name){1}) ||...
%             ~isfield(data{k}, indep_name)
%         error('The selected files do not match.')
%     end
% end

for k = 1:length(dep_vars)
    data_variable = dep_vars{k};

    for q = 1:length(filenames)
        % Check that the data variable exists (compute it if necessary).
        [data{q}, data_variable] = checkDataVar(data{q}, data_variable);
        dep_vals(q, :) = data{q}.(data_variable);
        dep_rels = data{q}.rels.(data_variable);
 
        if ~isfield(data{q}, data_variable)
            error(['Could not find data variable ''', data_variable, ''' in ',...
                'at least one of the files.'])
        elseif isempty(dep_rels)
                error(['Independent (sweep) variables for data variable ''',...
                      strrep(data_variable, '_', ' '),...
                      ''' are not specified.'])
        end
    end
    avg_data = mean(dep_vals);

%     processed_data_var = ['Average_', data_variable];
%     data2plot = data{q};
%     data2plot.(processed_data_var) = avg_data;
%     data2plot.units.(processed_data_var) = data{1}.units.(data_variable);
%     data2plot.rels.(processed_data_var) = {dep_rels{1}};
%     data2plot.dep{length(data{1}.dep) + 1} = processed_data_var;
%     data2plot.plotting.(processed_data_var).plot_title =...
%             {[strrep(filenames{1}, '_', '\_'), ' - ',...
%               strrep(filenames{end}, '_', '\_')],...
%              ['[', data{1}.Timestamp, ' - ',...
%               data{end}.Timestamp, ']']};
    figure
    t = linspace(0,5,6);
    p = plot(20*t, avg_data);
    xlabel('Idle Gate Time (us)', 'FontSize',14)   
    ylabel('I (V)','FontSize',14);
    title({[strrep(filenames{1}, '_', '\_'), ' - ',...
              strrep(filenames{end}, '_', '\_')],...
             ['[', data{1}.Timestamp, ' - ',...
              data{end}.Timestamp, ']']})
    
    p.Marker = '.'; 
    p.MarkerSize = 12;
    p.MarkerFaceColor = [0, 0.447 , 0.741];
    grid on
    axis tight
    [~, filename, ~] = fileparts(filenames{end});
    saveMeasData(data2plot, [filename, '_avg'])

    end
end