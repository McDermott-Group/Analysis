function plotPrepInfidelityEstimation
%plotPrepInfidelityEstimation  Plot the estimated upper bound for
%the preparation fidelity.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

probabilities = {'Pi_Pulse_Switching_Probability',...
                 'No_Pulse_Switching_Probability',...
                 'Dark_Switching_Probability'};
% Check that the data variable exists.
for k=1:3
    [data, probabilities{k}] = checkDataVar(data, probabilities{k});
    dep_rels = data.rels.(probabilities{k});
    if isempty(dep_rels)
        error(['Independent (sweep) variables for data variable ''',...
              strrep(data_variable, '_', ' '), ''' are not specified.'])
    end
end

p_pi = data.Pi_Pulse_Switching_Probability;
p_no = data.No_Pulse_Switching_Probability;
p_dark = data.Dark_Switching_Probability;

infidelity = (p_no - p_dark) ./ (p_pi - 2 * p_no);

processed_data_var = 'Preparation_Infidelity';
data.(processed_data_var) = infidelity;
data.units.(processed_data_var) = '';
data.rels.(processed_data_var) = dep_rels;
data.dep{length(data.dep) + 1} = processed_data_var;
if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    full_name = 'Ground State Preparation Infidelity';
else
    full_name = 'Excited State Preparation Infidelity';
end
data.plotting.(processed_data_var).full_name = full_name;

plotDataVar(data, processed_data_var);
end