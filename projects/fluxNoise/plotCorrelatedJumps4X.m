date = '02-27-20';
start_file = 0;
end_file = 829;

pathname = ['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Q3Q4Corr\General\' date '\QP_correlation100\MATLABData'];
% pathname = '/Volumes/smb/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q1Q2Corr/General/02-18-20/Charge_resetting/MATLABData';

n_jumps = zeros(1,end_file-start_file+1);
o_avg = zeros(4, end_file-start_file+1);

j = 1;
for i = start_file:end_file
    
    filename = ['QP_correlation100_',num2str(i,'%03d'),'.mat'];
    % Read the data file, convert the variable names, and specify the units.
    try
        data = cell(0);
        for k = 1:n
            file = fullfile(pathname, filename);
            data = loadMeasurementData(file);
        end
    catch
        error('The selected files contain unmatched data.')
    end
    
    if mod(j,100) == 0
        j
    end
    
    o1 = data.Single_Shot_Occupations_SB1(:);
    o2 = data.Single_Shot_Occupations_SB2(:);
    o3 = data.Single_Shot_Occupations_SB3(:);
    o4 = data.Single_Shot_Occupations_SB4(:);
    n_jumps(j) = sum(o1&o2&o3&o4);
    o_avg(:,j) = mean([o1 o2 o3 o4]);
    
    j = j + 1;

end

figure;
plot(n_jumps)
