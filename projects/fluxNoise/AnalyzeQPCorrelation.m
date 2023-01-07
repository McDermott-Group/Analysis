
dataType = 'QP_Correlation100';
ext = strcat(['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\02-17-20\',dataType,'\MATLABData\',dataType,'_']);
min_file_index = 0;
max_file_index = 980;

n_files = max_file_index - min_file_index + 1;
acf = zeros(21,1);
N_a_array = zeros(1,n_files*10);
N_b_array = zeros(1,n_files*10);
N_corr_array = zeros(1,n_files*10);
L = 10000;
N = n_files*10*L;
% a_all = zeros(1, N);
% b_all = zeros(1, N);
% ab_all = zeros(1, N);
for j = min_file_index:max_file_index
    ldata = loadMeasurementData([ext,num2str(j,'%03d'),'.mat']);
    a = ldata.Single_Shot_Occupations_SB1';
    b = ldata.Single_Shot_Occupations_SB2';
    ab = a(:) & b(:);
    L = size(a,1);
    
    i = j - min_file_index;
    
%     a_all(1+i*10*L:(i+1)*10*L) = a(:);
%     b_all(1+i*10*L:(i+1)*10*L) = b(:);
%     ab_all(1+i*10*L:(i+1)*10*L) = ab(:);
    
    N_a_array(10*i+1:10*i+10) = sum(a);
    N_b_array(10*i+1:10*i+10) = sum(b);
    N_corr_array(10*i+1:10*i+10) = sum(a&b);
    
    acf = acf + autocorrelate(ab)/n_files;
end
N = n_files*10*L;
N_a = sum(N_a_array);
N_b = sum(N_b_array);
N_corr = sum(N_corr_array);

[N_a/N, N_b/N]
randRate = N_a/N * N_b/N
corrRate = N_corr/N
figure; hold on; plot(0:20,acf); line([0 20],[corrRate^2 corrRate^2])
xlabel('Measurement Steps'); ylabel('Autocorrelation');

% high pass filter to remove LF oscillations in bg level
N_a_array = mean(N_a_array) - movmean(N_a_array, 100) + N_a_array;
N_b_array = mean(N_b_array) - movmean(N_b_array, 100) + N_b_array;
N_corr_array = mean(N_corr_array) - movmean(N_corr_array, 100) + N_corr_array;

figure; plot([N_a_array', N_b_array', N_corr_array']/L)

% figure; plot([movmean(a',10), movmean(b',10), movmean(ab',10)]);


function [acf] = autocorrelate(x)

acf = zeros(21,1);
for j = 1:21
    i = j-1;
    acf(j) = sum(x(1:end-i).*x(1+i:end))/(length(x)-i);
end

end