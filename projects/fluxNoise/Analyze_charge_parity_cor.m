
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

totFiles = maxFileIndex - minFileIndex + 1;
NumP = reps * trials;

cutoff = 200;
sum11ee = zeros(cutoff,1);
sum11eo = zeros(cutoff,1);
sum10ee = zeros(cutoff,1);
sum10eo = zeros(cutoff,1);
sum01ee = zeros(cutoff,1);
sum01eo = zeros(cutoff,1);
t = (1:cutoff)*100e-6;
nFiles = maxFileIndex - minFileIndex;
for i = minFileIndex:maxFileIndex
    ldata = load(strcat([ext,num2str(i,'%03d'),'.mat']));
    M = reshape(eval(strcat(['ldata.',dataType,'.Data.Single_Shot_Occupations_1']))',[reps*trials,1]);
    P = reshape(eval(strcat(['ldata.',dataType,'.Data.Single_Shot_Occupations_2']))',[reps*trials,1]);
    P = (P==M);
    
    for j = 1:cutoff
        dP = P(1:end-j)==P(j+1:end);
        Ma = M(1:end-j);  Mb = M(j+1:end);
        M0a = (Ma == 0);  M0b = (Mb == 0);
        M1a = (Ma == 1);  M1b = (Mb == 1);
        sum11ee(j) = sum11ee(j) + sum((M1a==1) + (M1b==1) + (dP==1) == 3)/cutoff/nFiles;
        sum11eo(j) = sum11eo(j) + sum((M1a==1) + (M1b==1) + (dP==0) == 3)/cutoff/nFiles;
        sum10ee(j) = sum10ee(j) + sum((M1a==1) + (M1b==0) + (dP==1) == 3)/cutoff/nFiles;
        sum10eo(j) = sum10eo(j) + sum((M1a==1) + (M1b==0) + (dP==0) == 3)/cutoff/nFiles;
        sum01ee(j) = sum01ee(j) + sum((M1a==0) + (M1b==1) + (dP==1) == 3)/cutoff/nFiles;
        sum01eo(j) = sum01eo(j) + sum((M1a==0) + (M1b==1) + (dP==0) == 3)/cutoff/nFiles;
    end
end

figure; plot(t,sum11ee,t,sum11eo,t,sum10ee,t,sum10eo,t,sum01ee,t,sum01eo); grid on;