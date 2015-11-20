function plotHistEqMeasData(data_variable)
%plotHistEqMeasData(DATA_VARIABLE)  Plot a histogram-equalized data.
%DATA_VARIABLE should be a name of the data variable to plot.

if ~exist('data_variable', 'var')
    error('No dependent data variable name is given as an input argument.')
end

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);
    
if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

dep_vals = (dep_vals - min(dep_vals(:))) ./...
           (max(dep_vals(:)) - min(dep_vals(:)));

processed_data_var = ['HistEqual_', data_variable];
data.(processed_data_var) = histeq(dep_vals);
data.units.(processed_data_var) = data.units.(data_variable);
data.rels.(processed_data_var) = data.rels.(data_variable);
data.dep{length(data.dep) + 1} = processed_data_var;
data.plotting.(processed_data_var).full_name =...
    ['Histogram-Equalized ', strrep(data_variable, '_', ' ')];

plotDataVar(data, processed_data_var);