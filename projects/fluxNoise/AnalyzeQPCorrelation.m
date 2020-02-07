
dataType = 'QP_Correlation464';
ext = strcat(['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\01-27-20\',dataType,'\MATLABData\',dataType,'_']);
min_file_index = 0;
max_file_index = 483;

n_files = max_file_index - min_file_index + 1;
acf = zeros(21,1);
N_a_array = zeros(1,n_files*10);
N_b_array = zeros(1,n_files*10);
N_corr_array = zeros(1,n_files*10);
for i = min_file_index:max_file_index
    ldata = loadMeasurementData([ext,num2str(i,'%03d'),'.mat']);
    a = ldata.Single_Shot_Occupations_a';
    b = ldata.Single_Shot_Occupations_b';
    
    N_a_array(10*i+1:10*i+10) = sum(a);
    N_b_array(10*i+1:10*i+10) = sum(b);
    N_corr_array(10*i+1:10*i+10) = sum(a&b);
    
    ab = a & b;
    ab = reshape(ab,10*10000,1);
    acf = acf + autocorrelate(ab)/n_files;
end
N = n_files*10*10000;
N_a = sum(N_a_array);
N_b = sum(N_b_array);
N_corr = sum(N_corr_array);

[N_a/N, N_b/N]
randRate = N_a/N * N_b/N
corrRate = N_corr/N
figure; hold on; plot(0:20,acf); line([0 20],[corrRate^2 corrRate^2])
xlabel('Measurement Steps'); ylabel('Autocorrelation');

figure; plot([N_a_array', N_b_array', N_corr_array']/10000)

function [acf] = autocorrelate(x)

acf = zeros(21,1);
for i = 1:21
    acf(i) = sum(x(1:end-i).*x(1+i:end))/(length(x)-i);
end

end