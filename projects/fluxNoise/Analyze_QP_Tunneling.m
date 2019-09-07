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
    
    [data_even, data_odd] = chop_data(data);
    
    [seg_psd_chopped] = crosspsd(data_even,data_odd,Fs/2);
    
    fileApsd_chopped = fileApsd_chopped + seg_psd_chopped(1:(NumP/4+1));
end
fileApsd_chopped = fileApsd_chopped * transfer / totFiles;

[windowApsd_chopped] = window_averaging(fileApsd_chopped);

%Fit to find QP tunneling rate:
% x = psd_freq_chopped(2:end);
% y = abs(windowApsd_chopped(2:end));
% [x, y] = prepareCurveData(x, y);
% [fr,~] = fitLorenzian(x,y);

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

function [psd] = crosspsd(z1,z2,Fs)
%Cross PSD
NumP = length(z1);
fft_seq1 = fft(z1);
fft_seq2 = conj(fft(z2));
psd = (1/(Fs*NumP))*(fft_seq1.*fft_seq2);
psd(2:(NumP/2)) = 2*psd(2:(NumP/2));
end

function [apsd] = window_averaging(psd)
% apsd = psd;
%Averaging with f window
apsd = zeros(size(psd));
for i = 1:size(psd,1)
    filter_fl = max(1,round(i - i/4));
    filter_fh = min(size(psd,1),round(i + i/4));
    apsd(i) = mean(psd(filter_fl:filter_fh));
end
end

function [seg_even, seg_odd] = chop_data(seg)
%Partition data into even and odd for CPSD on single data set
seg_even = seg(1:2:end);
seg_odd = seg(2:2:end);
end

function [fitresult, gof] = fitLorenzian(x, y)
% Set up fittype and options.
ft = fittype( 'a/(4*pi^2*x.^2+gamma^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [380 1056];

% Fit model to data.
[fitresult, gof] = fit( x, y, ft, opts );
end