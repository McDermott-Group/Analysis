path = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q4\General\03-06-20\Excess_1_state\MATLABData';

fs = 1/100e-6;
filesE1 = 7:106;
filesQP = 32:531;

% path = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q4\General\03-08-20\Excess_1_state\MATLABData';
% 
% fs = 1/100e-6;
% % filesE1 = [1:18 20:136 138:216 218:351 353:398 400:441 443:464 466:518];
% % filesE1 = [1341:1903];
% filesE1 = [0:436];
% filesQP = 32:531;

e1 = [];
i = 0;
for f = filesE1
    i = i + 1;
    try
        data = noiselib.load_file(path, ['Excess_1_state_',num2str(f,'%03d'),'.mat']);
        e1 = [e1 mean(data.Single_Shot_Occupations)];
    catch
        e1 = [e1 nan];
    end
end

figure; plot(e1)
% es = noiselib.unwrap_voltage_to_charge(data.Offset_Voltage, 0.2636, 1/0.2636); figure; hold on; plot(es); plot(30*e1); plot(es); plot(data.Offset_Voltage_R2)


highE1 = e1 >  0.1;
lowE1  = e1 <= 0.059;

highFileIndicies = filesQP(repelem(highE1,5));
lowFileIndicies  = filesQP(repelem(lowE1,5));

Plot_QP_Tunneling('DR1 - 2019-12-17\', 'CorrFar\', 'Q4\', '03-06-20\', highFileIndicies);
Plot_QP_Tunneling('DR1 - 2019-12-17\', 'CorrFar\', 'Q4\', '03-06-20\', lowFileIndicies);