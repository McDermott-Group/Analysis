function plotMeanNorm2DMeasData(data_variable, normalization_direction)
%plotMeanNorm2DMeasData(DATA_VARIABLE, NORMALIZATION_DIRECTION) Plot
%a line-by-line mean normalized 2D data.
%   plotMeanNorm2DMeasData(DATA_VARIABLE, NORMALIZATION_DIRECTION) plots
%   line-by-line rescaled data for vairable DATA_VARIABLE along
%   NORMALIZATION_DIRECTION. NORMALIZATION_DIRECTION should be either
%   'along_x' or 'along_y'.

if exist('normalization_direction', 'var') &&...
        ~isempty(strfind(normalization_direction, 'x'))
    normalization_direction = 'along_x';
else
    normalization_direction = 'along_y';
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

% Plot 2D data.
if length(dep_rels) == 2 
    if strcmp(normalization_direction, 'along_y')
        dep_vals = dep_vals ./ (mean(dep_vals, 2) * ones(1, size(dep_vals, 2)));
    elseif strcmp(normalization_direction, 'along_x')
        dep_vals = dep_vals ./ (ones(size(dep_vals, 1), 1) * mean(dep_vals));   
    end
    processed_data_var = ['MeanNorm_', data_variable];
    data.(processed_data_var) = dep_vals;
    data.units.(processed_data_var) = '';
    data.rels.(processed_data_var) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = processed_data_var;
    data.plotting.(processed_data_var).full_name =...
        ['Line-by-Line Mean-Normalized ', strrep(data_variable, '_', ' ')];
    data.plotting.(processed_data_var).extra_filename = ['_', normalization_direction];

    plotDataVar(data, processed_data_var);
end