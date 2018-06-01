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

ft = cell(length(data),1);
% nwashers = zeros(length(data),1);
qc = zeros(length(data),1);
qi = zeros(length(data),1);
startPoint = [100000,10000,0.95,5.53446,0,-45];
eqn = fittype('m*x + b + ( 1 - (Q0/Qc-2i*Q0*df/f0)/(1+2i*Q0*(f-f0)/f0) )*exp(1i*tau*2*pi*f)'); %Geerlings et al, Improving the quality factor of microwave compact resonators by optimizing their geometrical parameters
% eqn = fittype('m*x + b + 1/abs(1 + Qi/Qc*exp(j*phi)/(2*j*Qi*(x-f0)/f0))'); %Megrant et al, Planar superconducting resonators with internal quality factors above one million
figure(66)
set(gca,'fontsize', 18)
title({'5.5GHz Resonator','Output: 6 washers','Input: T-ed to make hanging-style'});
xlabel('Frequency [GHz]');
ylabel('Amplitude [Arb Units]');
for k = 1:length(data)
%     try
%         nwashers(k) = data{k}.Cavity_Input_Washers;
%     catch
%         nwashers(k) = 0;
%     end
    f = data{k}.RF_Frequency;
    s_dB = data{k}.S43;
    s_dB = s_dB - mean(s_dB);
    s = 10.^(s_dB/20);
    p = pi/180.*data{k}.S43_Phase;
    p = p - mean(p);
    i = s.*cos(p); 
    q = 1j*s.*sin(p);
    ss=(i+q);
    
    figure(22); plot(f, s_dB);
    figure(33); plot(f, 180/pi.*p);
    figure(44); plot(i, -1j*q)
    axis equal;
    figure(55); plot(real(1./ss),imag(1./ss))
%     ft{k} = fit(f.', s.', eqn, ...
%                 'StartPoint', startPoint, ...
%                 'Lower',[1,1,-Inf,5.5344,-Inf,-90], ...
%                 'Upper',[Inf,Inf,Inf,5.5345,Inf,90]);
    figure(66);
%     hold on;
    plot(f, s);
%     hold on;
%     plot(f, ft{k}(f), 'DisplayName', strcat('Fit: Q_c = ',num2str(ft{k}.Qc)));
%     figure(44);
%     hold on;
%     plot(ft{k}(f).'.*cos(p), ft{k}(f).'.*sin(p))
    
    opts = optimoptions(@lsqcurvefit,'Display','off');
%     objfcn = @(v,x)v(1)*x + v(2) + 1./(1 + v(3)/v(4)*exp(j*v(5))./(2*j*v(3)*(x-v(6))/v(6)));
% 'm*x + b + ( 1 - (Q0/Qc-2i*Q0*df/f0)/(1+2i*Q0*(f-f0)/f0) )*exp(1i*tau*2*pi*f'
    objfcn = @(v,x)[real(v(6)*x + v(2) + ( 1 - (v(1)/v(3)-2j*v(1)*v(4)/v(5))./(1+2j*v(1).*(x-v(5))/v(5)) ).*exp(1j*v(6)*2*pi*x)), ...
                    imag(v(6)*x + v(2) + ( 1 - (v(1)/v(3)-2j*v(1)*v(4)/v(5))./(1+2j*v(1).*(x-v(5))/v(5)) ).*exp(1j*v(6)*2*pi*x))];
%     objfcn = @(v,x)[real(v(1)*x + v(2) + 1./(1 + v(3)/v(4)*exp(j*v(5))./(2*j*v(3)*(x-v(6))/v(6)))), ... 
%                     imag(j*v(7)*x + j*v(8) + 1./(1 + j*v(9)/j/v(10)*exp(-v(11))./(2*j*v(9)*(x-j*v(12))/v(12))))];
%     v0 = [0; 0.95; 100000; 50000; -45; 5.53446; 1e-9; 1e-9; 1e-9; 1e-9; 1e-9; 1e-9];
    % Q0, b, Qc, df, f0, tau, m
    v0 = [50000; 0.6+0.25j; 50000; 0; 5.53446; 80e-9; 0];
    v0 = [6400; 2+5j; 10000; 0.00018; 6.7607; -80; 0];
%     [vestimated,resnorm] = lsqcurvefit(objfcn,v0,f.',ss.',[],[],opts)
    [vestimated,resnorm] = lsqcurvefit(objfcn,v0,f,[i,-1j*q],[],[],opts)
%     ft = objfcn(vestimated,f.');
    temp = objfcn(vestimated,f);
    ft = temp(:,1) + 1j*temp(:,2);
    
    figure(22); hold on;
    plot(f, 20*log10(abs(ft)))
    title('S21 in dB');
    figure(33); hold on;
    plot(f, 180/pi*atan(imag(ft)./real(ft)))
    title('S21 Phase');
    figure(44); hold on;
    plot(real(ft),imag(ft))
    title('IQ')
    figure(55); hold on;
    plot(real(1./ft),imag(1./ft))
    title('1/S21 IQ');
    figure(66); hold on;
    plot(f,abs(ft))
    title('S21');
end
legend('show');

% figure()
% plot(nwashers,qc)
% set(gca,'fontsize', 18)
% title('Extracted Q_c');
% xlabel('Number of washers');
% ylabel('Q_c');
% figure()
% plot(nwashers,qi)
% set(gca,'fontsize', 18)
% title('Extracted Q_i');
% xlabel('Number of washers');
% ylabel('Q_i');