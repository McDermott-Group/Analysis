
function data = processMeasurementData(data)
%processMeasurementData Rename some data fields, convert the units 
%if necessary, estimate errors and check that the data are properly
%specified.
%   DATA = processMeasurementData(DATA) renames some DATA fields,
%   converts the units, estimates errors when possible and check
%   the presence of some required fields. The fed DATA structure should
%   be taken from the output of the importMeasurementData function.

% Redefine the variable full names.
fields = fieldnames(data);
for k = 1:length(fieldnames(data))
    switch fields{k}
        case 'Init_Time'
            data = renameVariable(data, fields{k},...
                'Initialization_Time');

        case 'Probability'
            data = renameVariable(data, fields{k},...
                'Switching_Probability');

        case 'Number_of_Repetitions'
            if ~isnumeric(data.(fields{k}))
                data = renameVariable(data, fields{k},...
                    'Benchmarking_Number_of_Repetitions');
            end

        case 'Reps'
            if isfield(data, 'Actual_Reps') ||...
                    isfield(data, 'Number_of_Repetitions')
                data = renameVariable(data, fields{k},...
                    'Requested_Number_of_Repetitions');
            else
                data = renameVariable(data, fields{k},...
                    'Number_of_Repetitions');
            end

        case 'Actual_Reps'
            data = renameVariable(data, fields{k},...
                'Number_of_Repetitions');

        case 'Runs'
            data = renameVariable(data, fields{k},...
                'Number_of_Runs');

        case 'Pa'
            data = renameVariable(data, fields{k},...
                'JPM_A_Switching_Probability');

        case 'Pb'
            data = renameVariable(data, fields{k},...
                'JPM_B_Switching_Probability');

        case 'Detection_Time_Diff'
            data = renameVariable(data, fields{k},...
                'Detection_Time_Difference');

        case 'Detection_Time_Diff_Std_Dev'
            data = renameVariable(data, fields{k},...
                'Detection_Time_Difference_Std_Dev');

        case 'Corr_Coef'
            data = renameVariable(data, fields{k},...
                'Correlation_Coefficient');
    end
end

% Convert units.
fields = fieldnames(data);
for k = 1:length(fields)
    if isfield(data.units, fields{k})
        switch data.units.(fields{k})
            case 'DACUnits'
                data.units.(fields{k}) = 'DAC Units';
            case 'ADCUnits'
                data.units.(fields{k}) = 'ADC Units';
            case 'PreAmpTimeCounts'
                data.units.(fields{k}) = 'Preamp Time Counts';
        end
    end

    switch fields{k}
        case 'RF_Frequency'
            if isfield(data.units, fields{k}) &&...
                    strcmp(data.units.(fields{k}), 'Hz')
                data.units.(fields{k}) = 'GHz';
                data.(fields{k}) = data.(fields{k}) / 1e9;
            end

        case 'Readout_Frequency'
            if isfield(data.units, fields{k}) &&...
                    strcmp(data.units.(fields{k}), 'Hz')
                data.units.(fields{k}) = 'GHz';
                data.(fields{k}) = data.(fields{k}) / 1e9;
            end

        case 'Qubit_Frequency'
            if isfield(data.units, fields{k}) &&...
                    strcmp(data.units.(fields{k}), 'Hz')
                data.units.(fields{k}) = 'GHz';
                data.(fields{k}) = data.(fields{k}) / 1e9;
            end 

        case 'RF1_Frequency'
            if isfield(data.units, fields{k}) &&...
                    strcmp(data.units.(fields{k}), 'Hz')
                data.units.(fields{k}) = 'GHz';
                data.(fields{k}) = data.(fields{k}) / 1e9;
            end

        case 'RF2_Frequency'
            if isfield(data.units, fields{k}) &&...
                    strcmp(data.units.(fields{k}), 'Hz')
                data.units.(fields{k}) = 'GHz';
                data.(fields{k}) = data.(fields{k}) / 1e9;
            end
    end
end

% Estimate errors.
N = 1;
if isfield(data, 'Number_of_Repetitions')
    if ~isfield(data, 'Time') || length(data.Time) ~= 4096 ||...
            any(data.Time(:)' ~= 0:2:8190)
        N = double(N * data.Number_of_Repetitions);
    end
end
if isfield(data, 'Number_of_Runs')
    N = double(N * data.Number_of_Runs);
end
if N > 1
    N = N - 1;
end
for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if ~isfield(data.distr, dep_name)
        continue
    end 
    if strcmp(data.distr.(dep_name), 'binomial')
        data.error.(dep_name) = sqrt(double(data.(dep_name) .*...
            (1 - data.(dep_name)) / N));
    elseif strcmp(data.distr.(dep_name), 'normal')
        std_name = [dep_name, '_Std_Dev'];
        if isfield(data, std_name)
            data.error.(dep_name) = data.(std_name) / sqrt(N);
        else
            pos = strfind(dep_name, '_');
            if ~isempty(pos)
                ch = dep_name(pos(end):end);
                std_name = [dep_name(1:pos(end)), 'Std_Dev', ch];
                if isfield(data, std_name)
                    data.error.(dep_name) = data.(std_name) / sqrt(N);
                end
            end
        end
    % The following case is kept for the historical reasons.
    elseif strcmp(data.distr.(dep_name), 'uknown')
        err_name = [dep_name, '_Error'];
        if isfield(data, err_name)
            data.error.(dep_name) = data.(err_name);
        end
    elseif strcmp(data.distr.(dep_name), 'unknown')
        err_name = [dep_name, '_Error'];
        if isfield(data, err_name)
            data.error.(dep_name) = data.(err_name);
        end
    elseif strcmp(data.distr.(dep_name), 'std')
        continue
    end
end

% Sanity checks.
if isfield(data, 'Filename')
    selected_file = data.Filename;
elseif isfield(data, 'Experiment_Name') && isfield(data, 'Timestamp')
    selected_file = [data.Experiment_Name, ' [', data.Timestamp, ']'];
else
    selected_file = '[NO INFORMATION AVAILABLE]';
end
if isempty(data)
    error(['File ', selected_file, ' does not contain any data.']);
elseif ~isfield(data, 'indep') || isempty(data.indep)
    error(['The independent (sweep) variables are not specified ',...
        'in ', selected_file, '.']);
elseif ~isfield(data, 'dep') || isempty(data.dep)
    error(['The dependent (data) variables are not specified in ',...
        selected_file, '.']);
elseif ~isfield(data, 'rels') || isempty(data.rels)
    error(['The relationships between the dependent (data) and ',...
        'independent (sweep) variables are not specified in ',...
        selected_file, '.']);
else
    for k = 1:length(data.indep)
        if isempty(data.(data.indep{k}))
            error(['File ', selected_file, ' does not specify',...
                ' the independent (sweep) variables.']);
        end
    end
    for k = 1:length(data.dep)
        if isempty(data.(data.dep{k}))
            error(['File ', selected_file,...
                ' does not contain any actual data.']);
        end
    end
    RelationshipsFlag = false; 
    for k = 1:length(data.dep)
        if ~isempty(data.(data.dep{k}))
            RelationshipsFlag = true;
        end
    end
    if ~RelationshipsFlag
        error(['The relationships between the dependent (data) and',...
            ' independent (sweep) variables are not specified in ',...
            selected_file, '.']);
    end
end

% Delete dependent variables that don't contain any meaningful data.
for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    if all(~isfinite(data.(dep_name)))
        data = rmfield(data, dep_name);
        data.dep{data_index} = [];
        if isfield(data, 'rels') && isfield(data.rels, dep_name)
            data.rels = rmfield(data.rels, dep_name);
        end
        if isfield(data, 'distr') && isfield(data.distr, dep_name)
            data.distr = rmfield(data.distr, dep_name);
        end
        if isfield(data, 'units') && isfield(data.units, dep_name)
            data.units = rmfield(data.units, dep_name);
        end
        if isfield(data, 'error') && isfield(data.error, dep_name)
            data.units = rmfield(data.units, dep_name);
        end
    end
end

data.dep = data.dep(~cellfun('isempty', data.dep));

end