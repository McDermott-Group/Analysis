function data = plotFFT
%plotFFT  Plot FFT of the time domain data.
%   DATA = plotFFT The function returns structure DATA containing
%   the normalized data.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    dep_rels = data.rels.(I_name);
    
    if isempty(dep_rels) || length(dep_rels) > 2 ||...
            ~isempty(strfind(I_name, 'I')) ||...
            ~isempty(strfind(I_name, '_Std_Dev'))
        continue
    else
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data, Q_name)
            continue
        end
        for k = 1:length(dep_rels)
            if ~strcmp(dep_rels{k}, data.rels.(Q_name){k})
                continue
            end
        end
    end

    dep_vals = data.(I_name) + 1i * data.(Q_name);
    if isempty(strfind(dep_rels{1}, 'Time'))
        dep_vals = dep_vals';
        dep_rels = {dep_rels{2}, dep_rels{1}};
    end
    time = dep_rels{1};
    dt = median(diff(data.(time)));
    len = floor(length(data.(time)) / 2);
    fft_freq = 'FFT_Frequency';
    dep_rels{1} = fft_freq;
    data.(fft_freq) = linspace(1, len, len) /...
            (length(data.(time)) * dt);
    data.units.(fft_freq) = ['1/', data.units.(time)];
    
    fft_name = 'FFT_Spectrum';
    abs_fft = abs(fft(dep_vals));
    
    if length(dep_rels) == 1
        abs_fft = abs_fft(2:len + 1);
    else
        abs_fft = abs_fft(2:len + 1, :);
    end

    data.(fft_name) = abs_fft;
    data.units.(fft_name) = data.units.(I_name);
    data.rels.(fft_name) = dep_rels;
    data.dep{length(data.dep) + 1} = fft_name;
    data.plotting.(fft_name).full_name = '|FFT Spectrum|';
    data.plotting.(fft_name).extra_filename = '_fft';
     
    plotDataVar(data, fft_name);

end
end