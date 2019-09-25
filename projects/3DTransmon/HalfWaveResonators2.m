% Chris Wilen
% Probably not generally useful.  Varied number of washers on 3d cavity
% pins to see how that effecs the Qc.

% Select files to plot.
[filenames, pathnames, status] = selectMeasurementDataFile;
if ~status
    return
end

if ~iscell(filenames)
    filenames = {filenames};
end
if ~iscell(pathnames)
    pathnames = {pathnames};
end

% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.');
end

f = cell(length(data));
nwashers = zeros(length(data),1);
qc = zeros(length(data),1);
qi = zeros(length(data),1);
startPoint = [10000,2200,0.5,0.2435,-1.8,5.509];
startPoint = [100000,22000,0.5,0.2435,0,5.534525];
eqn = fittype('m*x + b + abs(1 - ((1/Qc)-2*i*dw/w0) / ((1/Qi+1/Qc)+2*i*(x-w0)/w0))'); %Geerlings et al, Improving the quality factor of microwave compact resonators by optimizing their geometrical parameters
% eqn = fittype('m*x + b + 1/(1 + Qi/Qc*exp(i*phi)/(2*i*Qi*(f-f0)/f0))'); %Geerlings et al, Improving the quality factor of microwave compact resonators by optimizing their geometrical parameters
figure()
set(gca,'fontsize', 18)
title({'5.5GHz Resonator','Output: 6 washers','Input: T-ed to make hanging-style'});
xlabel('Frequency [GHz]');
ylabel('Amplitude [Arb Units]');
for k = 1:length(data)
    try
        nwashers(k) = data{k}.Cavity_Input_Washers;
    catch
        nwashers(k) = 0;
    end
    f{k} = fit(data{k}.RF_Frequency, 10.^(data{k}.S43/10), eqn, ...
                'StartPoint', startPoint, ...
                'Lower',[1,1,-Inf,-Inf,-Inf,5.507], ...
                'Upper',[Inf,Inf,Inf,Inf,Inf,5.511]); 
    hold on;
    plot(data{k}.RF_Frequency, 10.^(data{k}.S43/10), 'DisplayName', strcat('Washers = ',num2str(nwashers(k))));
    hold on;
    plot(data{k}.RF_Frequency, f{k}(data{k}.RF_Frequency), 'DisplayName', strcat('Fit: Q_c = ',num2str(f{k}.Qc)));
    qc(k) = f{k}.Qc;
    qi(k) = f{k}.Qi;
    startPoint = coeffvalues(f{k});
    startPoint(1) = 1.75*startPoint(1);
end
legend('show');

figure()
plot(nwashers,qc)
set(gca,'fontsize', 18)
title('Extracted Q_c');
xlabel('Number of washers');
ylabel('Q_c');
figure()
plot(nwashers,qi)
set(gca,'fontsize', 18)
title('Extracted Q_i');
xlabel('Number of washers');
ylabel('Q_i');