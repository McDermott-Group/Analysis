function plotFPARFFreqData
%plotFPARFFreqData  Show increase in switching probability due to the RF
%excitation. The selected datafile should contain a 2D array of
%probability values with Fast Pulse Amplitude and RF Frequency in place of
%the independent variables.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isempty(strfind(dep_name, '_Std_Dev')) ||...
        isempty(strfind(dep_name, 'Probability'))
        continue
    end

    dep_vals = data.(dep_name);
    dep_rels = data.rels.(dep_name);
    
    if isempty(dep_rels)
        continue
    end

    % Plot 2D data.
    if length(dep_rels) == 2
        if strcmp(dep_rels{2}, 'Fast_Pulse_Amplitude')
            dep_vals = dep_vals - mean(dep_vals(:, 1:5), 2) *...
                ones(1, size(dep_vals, 2));
        else
            dep_vals = dep_vals - ones(size(dep_vals, 1), 1) *...
                mean(dep_vals(1:5, :));
        end

        processed_data_var = ['RFInduced_', dep_name];
        data.(processed_data_var) = dep_vals;
        data.units.(processed_data_var) = '';
        data.rels.(processed_data_var) = data.rels.(dep_name);
        data.dep{length(data.dep) + 1} = processed_data_var;
        data.plotting.(processed_data_var).full_name =...
            ['RF-Induced Increase in ', strrep(dep_name, '_', ' ')];

        plotDataVar(data, processed_data_var);
    end
end