function [nPhotons, Qi0, Qc0] = quarterWaveFitsLab1(data)

warning('off','all');
pause on;
att = 120;
gain = 70;

% data = loadMeasurementData
f = data.RF_Frequency;
logMagData = data.S21;
[m, index1] = min(logMagData);
fres = f(index1)
span = (max(f)-min(f))*1e3
thru = logMagData(1);
power = -att+data.NA_Source_Power


header = {'Photons','f0','df0','Qi','dQi','Qc','dQc'};


options = optimset('MaxFunEval',200000,'MaxIter',200000); 





%        f0                  qi          qc          L      thru      alpha
start = [fres,               5e3,      5e3,      0,     thru,    0,      0];
lower = [fres-span/4e3,    2.0e3,      1.0e3,     -50,    thru-50,     -50,    -100];
upper = [fres+span/4e3,    1.0e8,      1.0e6,      50,    10,      50,     100];

% Handle for fitting function mazinQuarter
mazinModel = @mazinQuarter;

% Handle for phase fitting function fitPhase
% %%%%%%%%%%%%phaseModel = @fitPhase;

% run fminseachbnd to calculate F0, Qi,and Qc
[mazinEstimates] = fminsearchbnd(mazinModel,start,lower,upper,options);


% Calculate confidence bounds on fit parameters
[error] = mazinError(f,logMagData,mazinEstimates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      Save fit estimates to local variables                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fFit = mazinEstimates(1);       % f0 fit
errorF = error(1);              % f0 error
qiFit = mazinEstimates(2);      % internal Q
errorQi = error(2);             % internal Q error
qcFit = mazinEstimates(3);      % Coupling Q
errorQc = error(3);             % Coupling Q error
lFit = mazinEstimates(4);       % Assymetry parameter
thruFit = mazinEstimates(5);    % through power level
aFit = mazinEstimates(6);       % Frequency dependant loss coeff.
zRFit = mazinEstimates(7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%               t      t0
tauStart      = [ 30e-9, 180];
tauLower      = [     0,   0];
tauUpper      = [100e-9, 360];

if power > 0
    power = -1*power;
end
if power<0
    power = power;
end
% convert between dBm and watts
% 1 mW = 0 dBm
powerWatts = 10^(power/10)/1000;
%calculate total quality factor
totalQ = (1/qiFit + 1/qcFit)^(-1);

% Planck's constant
hPlanck = 6.626e-34; % J*s

% Calcuate number of photons
nPhotons = powerWatts*(totalQ^2/qcFit/hPlanck/(fFit*1e9)^2);
Qc0 = qcFit;
Qi0 = qiFit;

% Save fits from loop iterations to array
dataFitResults(:) = [nPhotons,fFit,errorF,qiFit,errorQi,qcFit,errorQc];
% end

% dump array contents to output dataset
quarterWaveFitsLab1 = dataset({dataFitResults,header{:}});

    function S21 = mazinFunction(params,f)
        % resonant frequency
        f0 = params(1);
        % internal Q
        Qi = params(2);
        % coupling Q
        Qc = params(3);
        % Asymetry parameter
        L=params(4);
        % total Q
        Q=1/(1/Qc + 1/Qi);
        % min value past resonator
        Smin=Qc/(Qi+Qc);
        % reduced frequency
        dx=(f-f0)/f0;
        % S21 equation as defined by Mazin
        S21=(Smin + 2*1i*Q*dx)./(1 + 2*1i*Q*dx+1i*L);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert Mazin's S21 equation into LogMag    %%%
%%% while also accouting for loss in the        %%%
%%% cable                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s21LogMag = mazinFunctionLogMag(params,f)
        % Through power level (dB)
        through = params(5);
        % First frequency in file
        ff = f(1);
        % loss as a function of frequency (dB/Hz)
        alpha = params(6);
        % Calculate Mazin S21 signal
        signal= mazinFunction(params,f);
        % Apply corrections to Mazin S21
        s21LogMag = (alpha*(ff-f))+10*log10(signal.*conj(signal)) + through;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fitting equations for fminsearchbnd to extract  %%%
%%% f0.Qi,Qc                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [sse,error1,FittedCurve,signal] = mazinQuarter(params)
        % Through power level of resonator (dB)
        through = params(5);
        % first frequency in resonator file
        ff = f(1);
        % frequency dependant loss parameter (dB/Hz)
        alpha = params(6);
        % calculate Mazin S21 equation
        signal= mazinFunction(params,f);
        % alter signal to deal with frequency dependant loss and through level
        signal = (alpha*(ff-f))+ 10*log10(signal.*conj(signal)) + through;
        % save calculated signal to fitted curve
        FittedCurve = signal;
        % Compute the error vector (this is what fminsearchbnd minimizes)
        ErrorVector = FittedCurve - logMagData;
        % error1 is returned to fminseachbnd
        error1=ErrorVector;
        % sum of square error is returned to fminseachbnd
        sse = sum(ErrorVector.^ 2);
        % uncomment the following lines to watch fit in real time
        figure(14); clf;
        plot(f,logMagData,'b-','markersize',.5);    hold all;
        plot (f, signal, 'r-.','linewidth',1);
        drawnow; grid on; xlabel('RF Frequency (GHz)'); ylabel('S21 (dB)');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     The following code chunk calculates confidence bounds           %%%
%%%     for fit parameters                                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [paramError] = mazinError(f,logMagData,mazinEstimates)
        [paroutb,resid,J,Sigma] = nlinfit(f,logMagData,@mazinFunctionLogMag,mazinEstimates);
        ci = nlparci(paroutb,resid,'jacobian',J);
        paramError=abs(ci(:,2)-ci(:,1));
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     The following code plots the normalized data with mazin S21     %%%
%%%     in the complex plane so that fits can be analyzed               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotComplexFit(f0,qi,qc,L)
        % calc S21min defined by Mazin
        s21min = qc/(qc+qi);
        % calc total Q factor
        q = (qi*qc)/(qi+qc);
        %calc reduced frequency
        dx = (f-f0)/f0;
        %calc mazin S21 equation
        s21Eq = (s21min + 2*1i*q*dx)./(1 + 2*1i*q*dx+1i*L);
        %new figure
        figure(34);clf;
        % plot normalized data
        plot(real(normalData),imag(normalData),'bo');
        hold on;
        % plot fitted data
        plot(real(s21Eq),imag(s21Eq),'-r','linewidth', 3);
        % Wait for user input
        pause;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     The following code plots the logmag fit with params             %%%
%%%     calc'd from fminsearchbnd                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotLogMagFit(qi,qc,f0,l,through)
        
        s21min = qc/(qc+qi);
        q = (qi*qc)/(qi+qc);
        dx = (f-f0)/f0;
        s21Eq = (s21min + 2*1i*q*dx)./(1 + 2*1i*q*dx+1i*l);
        s21LogMag = 10*log10(s21Eq.*conj(s21Eq))+through;
        figure(35);clf;
        plot(f,logMagData,'b0');
        hold on;
        plot(f,s21LogMag,'-r','linewidth', 6);
        pause;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     The following code fits the phase parameters needed to          %%%
%%%     properly divide out the cable delay                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [sse,error1,FittedCurve] = fitPhase(params)
        % phi = f*tau + tau0
        tau = params(1);
        tau0 = params(2);
        signal = mazinFunction([fFit, qiFit,qcFit,lFit,thruFit,aFit,zRFit],f);
        signal = exp(-2*pi*1i*(f*tau + tau0)).*signal;
        signal = radtodeg(phase(signal));
        
        FittedCurve = signal;
        ErrorVector = FittedCurve - radtodeg(phase(complexData));
        error1=ErrorVector;
        sse = sum(ErrorVector.^ 2);
        % Uncomment the following lines to see phase fitting in real time
        %figure(15); clf;
        %plot(f,radtodeg(phase(complexData)),'bo','markersize',.5); hold all;
        %plot (f, signal, 'r--','linewidth',6);
        %drawnow;
    end



end