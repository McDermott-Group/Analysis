% % QP poisoning and excess 1 state
% date = '03-06-20\';
% filesE1 = 7:106;
% filesQP = 32:531;
% QPFilesPerE1 = 5;
% date = '03-11-20\'; % forgot to turn off during fill, too sparse to catch E1 jumps
% filesE1 = 1110:1204;
% filesQP = 0:473;
% QPFilesPerE1 = 5;
% date = '03-11-20\';
% filesE1 = 1212:1388;
% filesQP = 479:654;
% QPFilesPerE1 = 1;
% date = '03-11-20\';
% filesE1 = 1389:2381;
% filesQP = 655:1646;
% QPFilesPerE1 = 1;
date = '03-17-20\';
filesE1 = 0:917;
filesQP = 0:917;
QPFilesPerE1 = 1;

% % Charge offset and excess 1 state
% date = '03-08-20\';
% % filesE1 = [1:18 20:136 138:216 218:351 353:398 400:441 443:464 466:518];
% % filesE1 = [1341:1903];
% filesE1 = [0:436];
% date = '03-11-20\'; %cuy0632trn_parameters
% filesE1 = [110:1109];
% % filesE1 = [1110:1149];


fs = 1/100e-6;
% if we stopped the scan partyway through QP measurement
correctedNFilesE1 = fix(length(filesQP)/QPFilesPerE1);
filesE1 = filesE1(1:correctedNFilesE1); 
filesQP = filesQP(1:(correctedNFilesE1*QPFilesPerE1));
e1 = [];
acf = zeros(1,2*1000+1);
acf_baseline = zeros(1,2*1000+1);
n_high = 0; n_low = 0;
path = ['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1\General\' date 'Excess_1_state\MATLABData'];
for f = filesE1
    try
        data = noiselib.load_file(path, ['Excess_1_state_',num2str(f,'%03d'),'.mat']);
        o = data.Single_Shot_Occupations_SB1;
        e1_mean = mean(o);
        if e1_mean > 0.05
            n_high = n_high + 1;
            acf = (n_high-1)/n_high*acf + 1/n_high*noiselib.crosscorrelate(o, o, 1000);
        else
            n_low = n_low + 1;
            acf_baseline = (n_low-1)/n_low*acf_baseline + 1/n_low*noiselib.crosscorrelate(o, o, 1000);
        end
        e1 = [e1 e1_mean];
    catch
        e1 = [e1 nan];
    end
end

highThresh = 0.04;
lowThresh = 0.05;

figure; hold on; 
plot(e1); 
plot([0 length(e1)],[highThresh highThresh], 'k'); 
plot([0 length(e1)],[lowThresh lowThresh], 'k');
xlabel(['Files [' num2str(length(o(:))/fs) 'sec]']);
ylabel('Excess 1 State');

figure; hold on; plot(-1000:1000, acf_baseline); plot(-1000:1000, acf);
xlabel('Lags [100us]'); ylabel('Autocorrelation');
% es = noiselib.unwrap_voltage_to_charge(data.Offset_Voltage, 0.2636, 1/0.2636); figure; hold on; plot(es); plot(30*e1); plot(es); plot(data.Offset_Voltage_R2)


highE1 = e1 >  highThresh;
lowE1  = e1 <= lowThresh;

highFileIndicies = filesQP(repelem(highE1,QPFilesPerE1));
lowFileIndicies  = filesQP(repelem(lowE1,QPFilesPerE1));

Plot_QP_Tunneling('DR1 - 2019-12-17\', 'CorrFar\', 'Q1\', date, highFileIndicies);
Plot_QP_Tunneling('DR1 - 2019-12-17\', 'CorrFar\', 'Q1\', date, lowFileIndicies);