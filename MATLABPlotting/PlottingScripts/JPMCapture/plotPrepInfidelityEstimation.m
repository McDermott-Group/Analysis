function plotPrepInfidelityEstimation
%plotPrepInfidelityEstimation  Plot the estimated upper bound for
%the preparation fidelity.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

try
    probabilities = {'X_Gate_Switching_Probability',...
                     'I_Gate_Switching_Probability',...
                     'Dark_Switching_Probability'};
    % Check that the data variable exists.
    for k=1:length(probabilities)
        [data, probabilities{k}] = checkDataVar(data, probabilities{k});
        dep_rels = data.rels.(probabilities{k});
        if isempty(dep_rels)
            error(['Independent (sweep) variables for data variable ''',...
                  strrep(data_variable, '_', ' '), ''' are not specified.'])
        end
    end

    P_X = data.X_Gate_Switching_Probability;
    P_I = data.I_Gate_Switching_Probability;
    P_0 = data.Dark_Switching_Probability;

    % Check that the errors are given.
    error_flag = true;
    for k=1:length(probabilities)
        if ~isfield(data.error, probabilities{k}) 
            error_flag = false;
        end
    end

    if error_flag
        E_X = data.error.X_Gate_Switching_Probability;
        E_I = data.error.I_Gate_Switching_Probability;
        E_0 = data.error.Dark_Switching_Probability;
    end
catch
    probabilities = {'Pi_Pulse_Switching_Probability',...
                     'No_Pulse_Switching_Probability',...
                     'Dark_Switching_Probability'};
    % Check that the data variable exists.
    for k=1:length(probabilities)
        [data, probabilities{k}] = checkDataVar(data, probabilities{k});
        dep_rels = data.rels.(probabilities{k});
        if isempty(dep_rels)
            error(['Independent (sweep) variables for data variable ''',...
                  strrep(data_variable, '_', ' '), ''' are not specified.'])
        end
    end

    P_X = data.Pi_Pulse_Switching_Probability;
    P_I = data.No_Pulse_Switching_Probability;
    P_0 = data.Dark_Switching_Probability;

    % Check that the errors are given.
    error_flag = true;
    for k=1:length(probabilities)
        if ~isfield(data.error, probabilities{k}) 
            error_flag = false;
        end
    end

    if error_flag
        E_X = data.error.Pi_Pulse_Switching_Probability;
        E_I = data.error.No_Pulse_Switching_Probability;
        E_0 = data.error.Dark_Switching_Probability;
    end
end

if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    infidelity = (P_I - P_0) ./ (P_X - P_I);
    processed_data_var = 'Ground_State_Preparation_Infidelity';
    if error_flag
        E = infidelity .* (sqrt(E_I.^2 + E_0.^2) ./ (P_I - P_0) + ...
                           sqrt(E_X.^2 + E_I.^2) ./ (P_X - P_I));
    end
else
    infidelity = (P_X - P_0) ./ (P_X - P_I);
    processed_data_var = 'Excited_State_Preparation_Infidelity';
    if error_flag
        E = infidelity .* (sqrt(E_X.^2 + E_0.^2) ./ (P_X - P_0) + ...
                           sqrt(E_X.^2 + E_I.^2) ./ (P_X - P_I));
    end
end

E(imag(infidelity) ~= 0) = NaN;
infidelity(imag(infidelity) ~= 0) = NaN;
E(infidelity < 0) = NaN;
infidelity(infidelity < 0) = NaN;

sorted_infidelity = sort(infidelity);
sorted_infidelity(~isfinite(sorted_infidelity)) = [];
N = 20;
if length(sorted_infidelity) > N
    sorted_infidelity = sorted_infidelity(1:N);
end

percntiles = prctile(sorted_infidelity, [15 75]);
sorted_infidelity = sorted_infidelity(sorted_infidelity > percntiles(1) &...
                                      sorted_infidelity < percntiles(2));

bound = mean(sorted_infidelity);
bound_error = std(sorted_infidelity) / sqrt(N - 1);
if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    disp(['Ground State Preparation Infidelity = ', num2str(bound),...
        ' ± ', num2str(bound_error)]);
else
    disp(['Excited State Preparation Infidelity = ', num2str(bound),...
        ' ± ', num2str(bound_error)]);
end

data.(processed_data_var) = infidelity;
data.units.(processed_data_var) = '';
if error_flag
    data.error.(processed_data_var) = E;
end
data.rels.(processed_data_var) = dep_rels;
data.dep{length(data.dep) + 1} = processed_data_var;

plotDataVar(data, processed_data_var);

if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    infidelity = (1 - P_X) ./ (P_X - P_I);
    if error_flag
        E = infidelity .* (E_X ./ (1 - P_X) + ...
                sqrt(E_X.^2 + E_I.^2) ./ (P_X - P_I));
    end
    processed_data_var = 'Excited_State_Preparation_Infidelity';
else
    infidelity = (1 - P_I) ./ (P_I - P_X);
    if error_flag
        E = infidelity .* (E_I ./ (1 - P_I) + ...
                sqrt(E_I.^2 + E_X.^2) ./ (P_I - P_X));
    end
    processed_data_var = 'Ground_State_Preparation_Infidelity';
end

E(imag(infidelity) ~= 0) = NaN;
infidelity(imag(infidelity) ~= 0) = NaN;
E(infidelity < 0) = NaN;
infidelity(infidelity < 0) = NaN;

sorted_infidelity = sort(infidelity);
sorted_infidelity(~isfinite(sorted_infidelity)) = [];
N = 20;
if length(sorted_infidelity) > N
    sorted_infidelity = sorted_infidelity(1:N);
end

percntiles = prctile(sorted_infidelity, [15 75]);
sorted_infidelity = sorted_infidelity(sorted_infidelity > percntiles(1) &...
                                      sorted_infidelity < percntiles(2));

bound = mean(sorted_infidelity);
bound_error = std(sorted_infidelity) / sqrt(N - 1);
if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    disp(['Excited State Preparation Infidelity = ', num2str(bound),...
        ' ± ', num2str(bound_error)]);
else
    disp(['Groud State Preparation Infidelity = ', num2str(bound),...
        ' ± ', num2str(bound_error)]);
end

data.(processed_data_var) = infidelity;
data.units.(processed_data_var) = '';
if error_flag
    data.error.(processed_data_var) = E;
end
data.rels.(processed_data_var) = dep_rels;
data.dep{length(data.dep) + 1} = processed_data_var;

plotDataVar(data, processed_data_var);
end