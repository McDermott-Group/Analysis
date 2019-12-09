
minFileIndex = 3;
maxFileIndex = 650;

date = '10-15-19\';
samples = 'Circ1\';
qubit = 'Q3\';
dataType = 'Ramsey_XX';
CDdate = 'DR1 - 2019-09-24\';
ext = strcat(['Z:\mcdermott-group\data\fluxNoise\',CDdate,samples,qubit,'General\',date,dataType,'\MATLABData\',dataType,'_']);

nFiles = maxFileIndex - minFileIndex + 1;
df = zeros(nFiles,1);
o = zeros(nFiles,30);
o_fit = zeros(nFiles,30);
for i = 1:nFiles
    ldata = load(strcat([ext,num2str(i+minFileIndex,'%03d'),'.mat']));
    data = ldata.Ramsey_XX.Data;
    o(i,:) = data.Occupation;
    ft = fit(data.Idle_Gate_Time', data.Occupation', 'b + a*exp(-x/gamma)*cos(2*pi*df*x + phi)', ...
        'Start',[0.5 0.5 1.2e-3 0.4e4 0], 'Lower', [0.3 -1 1e-4 0.2e4 0]);
    df(i) = ft.df;
    o_fit(i,:) = ft(data.Idle_Gate_Time);
end

figure; imagesc(o);
figure; imagesc(o_fit);
figure; plot(df);