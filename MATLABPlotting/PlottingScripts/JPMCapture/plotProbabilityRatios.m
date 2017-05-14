function plotProbabilityRatios
%plotProbabilityRatios  Plot the switching probability ratios.

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
    ratio = P_X ./ P_0;
    processed_data_var = 'X_Gate_to_Dark_Ratio';
    if error_flag
        E = ratio .* (E_X ./ P_X + E_0 ./ P_0);
    end
else
    ratio = P_I ./ P_0;
    processed_data_var = 'I_Gate_to_Dark_Ratio';
    if error_flag
        E = ratio .* (E_I ./ P_I + E_0 ./ P_0);
    end
end

E(imag(ratio) ~= 0) = NaN;
ratio(imag(ratio) ~= 0) = NaN;
E(ratio < 0) = NaN;
ratio(ratio < 0) = NaN;

data.(processed_data_var) = ratio;
data.units.(processed_data_var) = '';
data.distr.(processed_data_var) = 'unknown';
if error_flag
    data.error.(processed_data_var) = E;
end
data.rels.(processed_data_var) = dep_rels;
data.dep{length(data.dep) + 1} = processed_data_var;

plotDataVar(data, processed_data_var);

if isfield(data, 'Driving_on_Dressed_One') &&...
        strcmp(data.Driving_on_Dressed_One, 'True')
    ratio = P_X ./ P_I;
    processed_data_var = 'X_Gate_to_I_Gate_Ratio';
    if error_flag
        E = ratio .* (E_X ./ P_X + E_I ./ P_I);
    end
else
    ratio = P_I ./ P_X;
    processed_data_var = 'I_Gate_to_X_Gate_Ratio';
    if error_flag
        E = ratio .* (E_I ./ P_I + E_X ./ P_X);
    end
end

E(imag(ratio) ~= 0) = NaN;
ratio(imag(ratio) ~= 0) = NaN;
E(ratio < 0) = NaN;
ratio(ratio < 0) = NaN;

data.(processed_data_var) = ratio;
data.units.(processed_data_var) = '';
data.distr.(processed_data_var) = 'unknown';
if error_flag
    data.error.(processed_data_var) = E;
end
data.rels.(processed_data_var) = dep_rels;
data.dep{length(data.dep) + 1} = processed_data_var;

plotDataVar(data, processed_data_var);
end