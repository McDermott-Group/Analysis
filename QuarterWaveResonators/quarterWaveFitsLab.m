function [quarterWaveFits] = quarterWaveFitsLab()

warning('off','all');
pause on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The goal of this code is to intake complex quarter wave resonator  %%%%
%%% data measured in transmission mode, convert it to logmag data, fit  %%%
%%% it using the equation derived by Mazin and then plot in the complex %%% 
%%% plane the fitted equation and the normalized complex data           %%%
%%%                                                                     %%%
%%%     File names should be in the following format...                 %%%
%%%     #.#### GHz Power -### dBm blahblahblah.s2p                       %%%
%%%                                                                     %%%
%%% This is because the center frequency and the measurement power      %%%
%%% are determined from the user supplied file name. Addtionally        %%%
%%% The measurement power is used in conjunction with fitted parameters %%%
%%% to calculate cavity photon occupation                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set header for data output   %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header = {'Photons','f0','df0','Qi','dQi','Qc','dQc'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% User defined Directory for data              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
complexDIR = uigetdir;  % Get the data directory

% Options for the number of times the given function for fminssearchbnd is
% to run
options = optimset('MaxFunEval',200000,'MaxIter',200000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load files from touchstone (s2p) %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
complexFiles = dir(fullfile(complexDIR,'*.s2p')); % Read all files with .s2p
 
% Define the # of rows for fitsReults to be number of files to read in
% This is where the fit results are stored
dataFitResults = zeros(length(complexFiles),7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The Fitting algorithm begins %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(complexFiles)
        
        % Read in the current file name
        name = complexFiles(i).name;
        
        % Set the resonant frequency from the first 6 file name characters
        % at this point the resonant frequency is in GHz
        fres = str2double(name(1:6));
        
        % Set the measurement power from the given information in file name
        power = str2double(name(18:21));
        
        % update the resonant frequency to be in units of Hz
        fres = fres*1e9; 
        
        % pointer to the current data file
        complexFileName = fullfile(complexDIR,complexFiles(i).name);
       
       
        % Read in ALL data contained in complexFileName
        complexData = dlmread(complexFileName,'',10,0);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
        % Frequency data is stored in the first column of complexData
        f = complexData(:,1);
        
        % Define which columns the real and complex portion of the data are
        % stored
        realDataColumn = 4;
        imagDataColumn = 5;
        
        % reduce compledData to just relevant data of measurement 
        complexData = complexData(:,realDataColumn) + 1i*complexData(:,imagDataColumn);
        
        % convert raw resonator data to log-magnitude for fitting
        logMagData = 10*log10(complexData.*conj(complexData));
       
        % convert raw resonator data to lin-magnitude for fitting
        linMagData = complexData.*conj(complexData);
        
        
        Thru = mean(logMagData(1:10));
        
        [~,j] = min(logMagData(:));
        
        fRes = f(j);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     the following chunk of code sets up and executes        %%%
        %%%     the fitting algorithm "fminsearchbnd.m" with            %%%
        %%%     initial guesses for fitting parameters in addition      %%%
        %%%     to upper and lower bounds on said fit parameteres       %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
             %   f0            qi     qc      L   thru  alpha
        start = [fRes,         2.0e4, 2.8e4,  0,   Thru,   0,   50];
        lower = [start(1)-5e6, 1.0e3, 1.1e3, -2,   Thru-10, -.5, -200];
        upper = [start(1)+5e6, 1.0e6, 1.0e6,  2,   Thru+10,  .5,  200];
        
        % Handle for fitting function mazinQuarter
        mazinModel = @mazinQuarter;
       
        % Handle for phase fitting function fitPhase
        phaseModel = @fitPhase;
       
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
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     The following fitting protocol extracts the added       %%%
        %%%     cable phase delay (t and t0) in the following form      %%%
        %%%
        %%%     measuredSignal = exp[2*pi*i*(f*t + phi_0)]*signal          %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %               t      t0
        tauStart      = [ 30e-9, 180];
        tauLower      = [     0,   0];
        tauUpper      = [100e-9, 360];
        
        % Fit t and t0 with phase model handle
        [tauEstimates] = fminsearchbnd(phaseModel,tauStart,tauLower,tauUpper,options);
        % save fit estimates to local variables
        timeDelay = tauEstimates(1);
        phi0 = tauEstimates(2);
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     In complex space, cable delay shows up as circles       %%%
        %%%     that are a function of the frequency. The following     %%%
        %%%     code is meant to extract the radius of the cable delay  %%%
        %%%     circles in phase space with a 2 parameter fit. One      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % We take just the first and last 1/10 of freq/complex data
        % to fit the radius
        [radiusFits] = polyfit([f(1:end/10);f(9*end/10:end)],[abs(complexData(1:end/10));abs(complexData(9*end/10:end))],1);
        
        % The radius may not be constant (due to loss) thus the two
        % parameter fit
        radius = f.*radiusFits(1)+radiusFits(2);
        
        % We now have everything we need to correctly model the
        % cable phase evolution in absence of the resonator
        phaseEvolution = radius.*exp(-2*pi*1i*(f.*timeDelay + phi0));
        
        % The normalization of the raw data simply involves dividing
        % the raw complex data by the fitted phase evolution 
        normalData = complexData./phaseEvolution;
        
        % Uncomment the following line to see the normalized data
        % with the complex fit
        % plotComplexFit(fFit,qiFit,qcFit,lFit)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Convert power input to cavity photon number             %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %insure that negative was correctly added into power
        if power > 0
            power = -1*power;
        end
        
        % convert between dBm and watts
        % 1 mW = 0 dBm
        powerWatts = 10^(power/10)/1000; 
        %calculate total quality factor
        totalQ = (1/qiFit + 1/qcFit)^(-1);
        % Planck's constant
        hPlanck = 6.626e-34; % J*s
        % Calcuate number of photons
        nPhotons = powerWatts*(totalQ^2/qcFit/hPlanck/fFit^2);
        
        % Save fits from loop iterations to array
        dataFitResults(i,:) = [nPhotons,fFit,errorF,qiFit,errorQi,qcFit,errorQc];
    end
    
    % dump array contents to output dataset
    quarterWaveFits = dataset({dataFitResults,header{:}});
        
function S21 = mazinFunction(params,f)
    % resonant frequency
    f0 = params(1);
    % internal Q
    Qi = params(2);
    % coupling Q
    Qc = params(3);
    % Assymetry parameter
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
%     figure(14); clf;
%     plot(f,logMagData,'b-','markersize',.5);    hold all;
%     plot (f, signal, 'r--','linewidth',6);
%     drawnow;
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
    plot(f,s21LogMag,'-r','linewidth', 3);
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