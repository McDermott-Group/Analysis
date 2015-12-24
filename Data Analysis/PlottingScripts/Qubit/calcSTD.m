function calcSTD
%calcSTD  Calculate standard deviation of the amplitude for a single short
%experiment data.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    Q_name = strrep(I_name, 'I', 'Q');
    
    if ~strcmp(I_name, 'Is') || ~isfield(data, Q_name)
        continue
    else
        Is = data.(I_name)(:);
        Qs = data.(Q_name)(:);
    end
    
    amplitude = sqrt((Is - mean(Is)).^2 + (Qs - mean(Qs)).^2);
    
    amplitude_std = std(amplitude);
    units = getUnits(data, I_name);
    
    disp(['Standard deviation of the amplitude = ',...
        num2str(amplitude_std), units, '.'])
end
end