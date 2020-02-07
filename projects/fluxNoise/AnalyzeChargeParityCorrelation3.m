% check QP tunneling rates 


%PARAMETERS
reps = 8192;
trials = 5;
minFileIndex = 0;
maxFileIndex = 585; 1113;
refreshTime = 100e-6;

date = '09-15-19\';
samples = 'Circ1\';
qubit = 'Q4\';
dataType = 'QP_Tunneling_PSD';
CDdate = 'DR1 - 2019-08-12\';
ext = strcat(['Z:\mcdermott-group\data\fluxNoise\',CDdate,samples,qubit,'General\',date,dataType,'\MATLABData\',dataType,'_']);

NumP = reps * trials;

cutoff = 200;
tot = zeros(cutoff,1);
t = (1:cutoff)*100e-6;
for i = minFileIndex:maxFileIndex
    ldata = load(strcat([ext,num2str(i,'%03d'),'.mat']));
    P = reshape(eval(strcat(['ldata.',dataType,'.Data.Single_Shot_Occupations']))',[reps*trials,1]);
    
    for j = 1:cutoff
        tot(j) = tot(j) + sum( P(1:end-j)==P(j+1:end) )/(NumP-j)/(maxFileIndex-minFileIndex);
    end
end

figure; plot(t,tot-0.5); grid on;