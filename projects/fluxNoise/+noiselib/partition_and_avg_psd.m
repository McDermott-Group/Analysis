function [avg_cpsd, cpsd_freq] = partition_and_avg_psd(data, Fs)
    % Splits the data into odd and even segments and takes the self cpsd.
    % data is input as a matrix, where each line (the second axis) will
    % give a cpsd that is then averaged along the first axis.
    N = size(data,2); % number of points per line
    n = size(data,1); % number of lines to average
    avg_cpsd = zeros(1,N/4+1);
    for i = 1:n
        [seg_even, seg_odd] = noiselib.partition_data(data(i,:));
        [seg_cpsd] = noiselib.crosspsd(seg_even,seg_odd,Fs/2);
        avg_cpsd = avg_cpsd + seg_cpsd(1:(N/4+1));
    end
    avg_cpsd = avg_cpsd/n;
    cpsd_freq = 0:Fs/N:(Fs/4);
end