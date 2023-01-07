% Copyright Chris Wilen 2018
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
    
    opts = optimoptions(@lsqcurvefit,'Display','off','MaxIterations',10000,'MaxFunEvals',10000, 'TolX', 1e-9,'TolFun',1e-9, 'StepTolerance',1e-9);%, 'PlotFcn',@optimplotresnorm);
%     objfcn = @(v,x)v(1)*x + v(2) + 1./(1 + v(3)/v(4)*exp(j*v(5))./(2*j*v(3)*(x-v(6))/v(6)));
% 'm*x + b + ( 1 - (Q0/Qc-2i*Q0*df/f0)/(1+2i*Q0*(f-f0)/f0) )*exp(1i*tau*2*pi*f+1i*phi0 )'
%     % Q0, Qc, df, f0, tau, b_re, m, b_im, phi
%     v0 = [370000; 450000; 0.05; 5.53446; -80e-9; 0; .2; 0; 0.5];
%     [vestimated,resnorm] = lsqcurvefit(objfcn,v0,f.',ss.',[],[],opts)
%     ft = objfcn(vestimated,f.');
    objfcn = @(v,x)[real(v(10)*( 1 - (v(1)/v(2)-2j*v(1)*v(3)/v(4)) ./ (1+2j*v(1).*(x-v(4))/v(4)) ).* exp(1j*(2*pi*v(5)*(x-v(4))+v(6))) ) + 1*v(7), ...
                    imag(v(10)*( 1 - (v(1)/v(2)-2j*v(1)*v(3)/v(4)) ./ (1+2j*v(1).*(x-v(4))/v(4)) ).* exp(1j*(2*pi*v(5)*(x-v(4))+v(6))) ) + 1*v(8)];
    % Q0, Qc, df, f0, tau, phi, b_re, b_im, m, scale factor
    v0 = [350000; 450000; 0.05; 5.53446; -80e-9; 0; mean(i); mean(q); 0.5; 1];
    % v0 = [2000; 3000; 0.05; 7.5466; -80e-9; 0; mean(i); mean(q); 0.5; 1];
    [vestimated,resnorm] = lsqcurvefit(objfcn,v0,f.',[i.',-1j*q.'],[],[],opts)
    1/(1/vestimated(1) - 1/vestimated(2))
    temp = objfcn(vestimated,f.');
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