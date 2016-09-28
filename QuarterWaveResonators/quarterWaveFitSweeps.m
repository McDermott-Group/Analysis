function [quarterWaveFits] = quarterWaveFitSweeps()

    warning('off','all');
    pause on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The goal of this code is to intake logmag data as taken via         %%%
%%% parameter (flux, power, etc) sweeps and fit the S21 function as     %%%
%%% derived by Mazin. This function requires the dependancy             %%%
%%% fminsearchbnd.m                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hPlanck = 6.626e-34;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     User data load in and independant variable determination        %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Load measurement data from sweep file
    data = loadMeasurementData;
    
    % Determine what is being swept
    sweptParamString = data.indep{1};
    
    % Get values for parameter that is being swept
    sweptParamArray = getfield(data,sweptParamString);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Set header for data output   %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    header = {sweptParamString,'f0','df0','Qi','dQi','Qc','dQc'};
    
    % Mutual inductance, H
    fluxMutual = 2.2e-12;
    
    % flux quantum
    phi0 = 2.068e-15;
    
    % Cavity Input Attenutation (dB)
    inputAttenuation = -100;
    
    % Set options for fitting routine
    options = optimset('MaxFunEval',200000,'MaxIter',200000);
    
    % set the frequency span for the data sets
    freqData = data.RF_Frequency*1e9;
    
    % Number of loop iterations
    loopEnd = length(data.S21(:,1));
    
    % Declare results bin loopEnd X length(header)
    dataFitResults = zeros(loopEnd,length(header));
    
    % Handle for fitting function mazinQuarter
    mazinModel = @mazinQuarter;
    
    % Loop over data sets and fit each logMag profile
    for i = 1:loopEnd
    
        logMagData = data.S21(i,:)';
        %figure(1000);clf; plot(freqData,logMagData); hold on;
        [~,fMinIndex] = min(logMagData);
        
        fRes = freqData(fMinIndex);
        thru = mean(logMagData(1:10));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     the following chunk of code sets up and executes        %%%
        %%%     the fitting algorithm "fminsearchbnd.m" with            %%%
        %%%     initial guesses for fitting parameters in addition      %%%
        %%%     to upper and lower bounds on said fit parameteres       %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %        f0            qi     qc      L   thru  alpha
        start = [fRes,         2.0e4, 2.8e4,  0,   thru,   0,   50];
        lower = [start(1)-5e6, 5.0e2, 1.1e3, -2,   thru-10, -.5, -200];
        upper = [start(1)+5e6, 1.0e6, 1.0e5,  2,   thru+10,  .5,  200];
        
        %freqHandle = figure('Visible','off');
        %setappdata(freqHandle,'freq',freq);
        
        %logMagHandle = figure('Visible','off');
        %setappdata(logMagHandle,'logMagData',logMagData);
        [mazinEstimates] = fminsearchbnd(mazinModel,start,lower,upper,options);
        
        % Calculate confidence bounds on fit parameters
        [error] = mazinError(freqData,logMagData,mazinEstimates);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%      Save fit estimates to local variables                  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fFit = mazinEstimates(1);       % f0 fit
        errorF = error(1);              % f0 error
        qiFit = mazinEstimates(2);      % internal Q
        errorQi = error(2);             % internal Q error
        qcFit = mazinEstimates(3);      % Coupling Q
        errorQc = error(3);             % Coupling Q error
        
        % Save fits from loop iterations to array
        dataFitResults(i,:) = [sweptParamArray(i),fFit,errorF,qiFit,errorQi,qcFit,errorQc];
        
    end
    
    plotData
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to determine what sweep independent    %%%
%%% variable is and plot internal Q                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotData
        if strcmp(sweptParamString,'Flux_Bias_Current_2') == 1 ||... 
                strcmp(sweptParamString,'Flux_Bias_Current_1') == 1
            % sweptParamArray = sweptParamArray / fluxBiasResistance;
            sweptParamArray = sweptParamArray* 1e-3 * fluxMutual;
            sweptParamArray = sweptParamArray / phi0;
            header{1} = 'Flux Bias, \Phi/\Phi_0';
            dataFitResults(:,1) = sweptParamArray; 
            % dump array contents to output dataset
            dataFitResults = sortrows(dataFitResults);
            quarterWaveFits = dataset({dataFitResults,header{:}});

            figure(4001); clf;
            errorbar(dataFitResults(:,1),dataFitResults(:,4),dataFitResults(:,5),'bx');
            ylabel('Internal Quality Factor');
            xlabel(header{1});
        elseif strcmp(sweptParamString,'NA_Source_Power') == 1
            sweptParamArray = sweptParamArray + inputAttenuation;
            %sweptParamArray = sweptParamArray / fluxBiasResistance;
            powerWatts = 10.^(sweptParamArray/10)/1000;
            fs = dataFitResults(:,2);
            Qis = dataFitResults(:,4);
            Qcs = dataFitResults(:,6);
            totalQ = (1./Qis + 1./Qcs).^(-1);
            nPhotons = powerWatts.*((totalQ.^2)./Qcs./hPlanck./fs.^2);
            sweptParamArray = nPhotons;
            header{1} = 'Cavity Photons $\bar{n}$';
            dataFitResults(:,1) = sweptParamArray; 
            % dump array contents to output dataset
            dataFitResults = sortrows(dataFitResults);
            quarterWaveFits = dataset({dataFitResults,header{:}});

            figure(4000); clf;

            errorbar(dataFitResults(:,1),dataFitResults(:,4),dataFitResults(:,5),'bx');
            ylabel('Internal Quality Factor');
            h = xlabel(header{1});
            set(h,'interpreter','Latex','FontSize',18);
            set(gca,'XScale','log')
        else
            % dump array contents to output dataset
            dataFitResults = sortrows(dataFitResults);
            quarterWaveFits = dataset({dataFitResults,header{:}});

            figure(3999); clf;
            errorbar(dataFitResults(:,1),dataFitResults(:,4),dataFitResults(:,5),'bx');
            ylabel('Internal Quality Factor');
            xlabel(header{1});
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% S21 function as defined in Ben Mazin's Thesis   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function S21 = mazinFunction(params,freqData)
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
        dx=(freqData-f0)/f0;
        % S21 equation as defined by Mazin
        S21=(Smin + 2*1i*Q*dx)./(1 + 2*1i*Q*dx+1i*L);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fitting equations for fminsearchbnd to extract  %%%
%%% f0.Qi,Qc                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [sse,error1,FittedCurve,signal] = mazinQuarter(params)
        % Through power level of resonator (dB)
        through = params(5);
        % first frequency in resonator file
        ff = freqData(1);
        % frequency dependant loss parameter (dB/Hz)
        alpha = params(6);
        % calculate Mazin S21 equation
        signal= mazinFunction(params,freqData);
        % alter signal to deal with frequency dependant loss and through level
        signal = (alpha*(ff-freqData))+ 10*log10(signal.*conj(signal)) + through;
        % save calculated signal to fitted curve
        FittedCurve = signal;
        % Compute the error vector (this is what fminsearchbnd minimizes)
        ErrorVector = FittedCurve - logMagData;
        % error1 is returned to fminseachbnd
        error1=ErrorVector;
        % sum of square error is returned to fminseachbnd
        sse = sum(ErrorVector.^ 2);
        % uncomment the following lines to watch fit in real time
        %figure(14); clf;
        %plot(freqData,logMagData,'b-','markersize',.5);    hold all;
        %plot (freqData, signal, 'r--','linewidth',6);
        %drawnow;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert Mazin's S21 equation into LogMag    %%%
%%% while also accouting for loss in the        %%%
%%% cable                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function s21LogMag = mazinFunctionLogMag(params,freqData)
        % Through power level (dB)
        through = params(5);
        % First frequency in file
        ff = freqData(1);
        % loss as a function of frequency (dB/Hz)
        alpha = params(6);
        % Calculate Mazin S21 signal
        signal= mazinFunction(params,freqData);
        % Apply corrections to Mazin S21
        s21LogMag = (alpha*(ff-freqData))+10*log10(signal.*conj(signal)) + through;
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


end









