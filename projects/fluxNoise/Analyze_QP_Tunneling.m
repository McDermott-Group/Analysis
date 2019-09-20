function [psd_freq_chopped,windowApsd_chopped] = Analyze_QP_Tunneling(samples, qubit, date, minFileIndex, maxFileIndex)
%Returns:
% windowApsd_charge_chopped: Charge self-CPSD
% windowApsd_flux_chopped: Flux self-CPSD
% windowAcpsd: Charge-flux CPSD
% psd_freq_chopped: Frequencies for self-CPSDs
% psd_freq: Frequencies for CPSD
% totFiles: Total files


%PARAMETERS
reps = 8192;
trials = 5;
% minFileIndex = 0;
% maxFileIndex = 83;
refreshTime = 100e-6;

% date = '08-22-19\';
% samples = 'Circ1\';
% qubit = 'Q2\';
dataType = 'QP_Tunneling_PSD';
% CDdate = 'DR1 - 2019-06-10\';
CDdate = 'DR1 - 2019-08-12\';
ext = strcat(['Z:\mcdermott-group\data\fluxNoise\',CDdate,samples,qubit,'General\',date,dataType,'\MATLABData\',dataType,'_']);

totFiles = maxFileIndex - minFileIndex + 1;
Fs = 1/refreshTime;%samples per second
NumP = reps * trials;
vis = 1;
transfer = 1/vis^2;%1.485e-5/0.0002957;
psd_freq_chopped = 0:Fs/NumP:(Fs/4);

fileApsd_chopped = zeros(reps*trials/4+1,1);
for i = minFileIndex:maxFileIndex
    ldata = load(strcat([ext,num2str(i,'%03d'),'.mat']));
    data = reshape(eval(strcat(['ldata.',dataType,'.Data.Single_Shot_Occupations']))',[reps*trials,1]);
    
    [data_even, data_odd] = noiselib.chop_data(data);
    
    [seg_psd_chopped] = noiselib.crosspsd(data_even,data_odd,Fs/2);
    
    fileApsd_chopped = fileApsd_chopped + seg_psd_chopped(1:(NumP/4+1));
end
fileApsd_chopped = fileApsd_chopped * transfer / totFiles;

[windowApsd_chopped] = noiselib.window_averaging(fileApsd_chopped);

%Fit to find QP tunneling rate:
% x = psd_freq_chopped(2:end);
% y = abs(windowApsd_chopped(2:end));
% [x, y] = prepareCurveData(x, y);
% [fr,~] = noiselib.fit_lorenzian(x,y);

% f = 0.01:0.01:5e3; gamma = 1056;lor = 3.8e2./(4*pi^2*f.^2+gamma^2);
% figure(12);hold on;plot(f,lor)

%Plot
figure(111);hold on
title('1/f Averaged PSD')
plot(psd_freq_chopped,abs(windowApsd_chopped), 'DisplayName', [samples(1:end-1),qubit(1:end-1)])
% plot(fr, x, y);
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Frequency (Hz)')
ylabel('S_\eta (\eta^2/Hz)')
grid on

end