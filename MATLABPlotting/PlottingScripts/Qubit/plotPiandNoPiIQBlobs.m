function [data] = plotPiandNoPiIQBlobs(GATE_I, GATE_X, PLOT_COMBINED)
%plotPiandNoPiIQBlobs(GATE_I, GATE_X, PLOT_COMBINED) plot IQ blobs from
%   many trials.
%
% Functions to plot two identical experiments where the difference is that
%   one of the datasets is from applying an idle gate and the other is from
%   applying a pi gate.
%
% Plots the dots and also returns stacked data in a struct.

if ~exist('GATE_I','var')
    % Select files to compute the difference.
    [filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the file corresponding to the Idle gate...',...
     'Select the file corresponding to the Pi gate...'});
    if ~status
        return
    end

    % Read the data files, convert the variable names, and specify the units.
    fileI = fullfile(pathnames{1}, filenames{1});
    fileX = fullfile(pathnames{2}, filenames{2});
    gateI = processMeasurementData(importMeasurementData(fileI));
    gateX = processMeasurementData(importMeasurementData(fileX));
else
    gateI = GATE_I;
    gateX = GATE_X;
end

if ~exist('PLOT_COMBINED','var')
    plotCombined = 0;
else
    plotCombined = PLOT_COMBINED;
end

%% Initialize and stack the data
data.gateI.Is = [];
data.gateI.Qs = [];
data.gateX.Is = [];
data.gateX.Qs = [];

for n = 1:length(gateI.Trial)
   data.gateI.Is = [data.gateI.Is gateI.Is(n,:)];
   data.gateI.Qs = [data.gateI.Qs gateI.Qs(n,:)];
   data.gateX.Is = [data.gateX.Is gateX.Is(n,:)];
   data.gateX.Qs = [data.gateX.Qs gateX.Qs(n,:)];
end

data.gateI.IQs = [data.gateI.Is' data.gateI.Qs'];
data.gateX.IQs = [data.gateX.Is' data.gateX.Qs'];

data.all.Is = [data.gateI.Is(1,:) data.gateX.Is(1,:)]';
data.all.Qs = [data.gateI.Qs(1,:) data.gateX.Qs(1,:)]';
data.all.IQs = [data.all.Is data.all.Qs];

%% Calculate mean IQ values via real data and via kmeans clustering
% mean values
data.gateI.I = mean(data.gateI.Is);
data.gateI.Q = mean(data.gateI.Qs);

data.gateX.I = mean(data.gateX.Is);
data.gateX.Q = mean(data.gateX.Qs);

% kmeans clustering
opts = statset('Display','final',...
               'UseParallel', 0,...
               'MaxIter', 10000);
[data.gateI.idx, data.gateI.C] = kmeans(data.gateI.IQs, 1,...
                                            'Replicates', 1,...
                                            'Options', opts);
[data.gateX.idx, data.gateX.C] = kmeans(data.gateX.IQs, 1,...
                                            'Replicates', 1,...
                                            'Options', opts);

if plotCombined == 1
    [data.all.idx, data.all.C] = kmeans(data.all.IQs, 2, 'Replicates', 1,...
                                                     'Options', opts);

    %Find bisecting line
    data.all.bslope = -((data.all.C(2,2) - data.all.C(1,2))/...
                       (data.all.C(2,1) - data.all.C(1,1)))^-1;
    data.all.byint = mean(data.all.C(:,2)) - data.all.bslope * mean(data.all.C(:,1));

    data.all.angle = atan(data.all.bslope);

    data.all.R = [cos(-data.all.angle) -sin(-data.all.angle);...
                  sin(-data.all.angle)  cos(-data.all.angle)];

    data.all.IQsTransformed = data.all.IQs * data.all.R;
end

          
data.trans.bslope = -((data.gateX.C(2) - data.gateI.C(2))/...
                   (data.gateX.C(1) - data.gateI.C(1)))^-1;
data.trans.byint = mean(data.gateX.C(2)) - data.trans.bslope * mean(data.gateX.C(1));

data.trans.angle = atan(data.trans.bslope);

data.trans.R = [cos(data.trans.angle) -sin(data.trans.angle);...
              sin(data.trans.angle)  cos(data.trans.angle)];

data.gateX.IQsTransformed = data.gateX.IQs * data.trans.R;
data.gateI.IQsTransformed = data.gateI.IQs * data.trans.R;


%% Do the plotting
% mean values
if plotCombined == 1
    figure;
    subplot(1,2,1)
    plot(data.gateI.Is, data.gateI.Qs,'b.',...
        'MarkerSize',2,'DisplayName','0')
    hold on
    plot(data.gateX.Is, data.gateX.Qs,'r.',...
        'MarkerSize',2,'DisplayName','1')
    xlabel('I (V)')
    ylabel('Q (V)')
    title({'Raw Quadrature Space Plot for ', ...
        [num2str(length(data.gateI.Is)) ' Shots']})
    axis square
    grid on
    set(gca,'FontSize',14)
    plot([data.gateI.I data.gateX.I],[data.gateI.Q data.gateX.Q],'g+',...
        'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');
    legend('show')
    hold off

    % kmeans clustering
    if data.all.idx(1) == 1
        blobAColor = 'b.';
        blobBColor = 'r.';
    else
        blobAColor = 'r.';
        blobBColor = 'b.';
    end
    subplot(1,2,2)
    plot(data.all.IQs(data.all.idx==1,1),data.all.IQs(data.all.idx==1,2),blobAColor,...
        'MarkerSize',2,'DisplayName','Blob A')
    hold on
    plot(data.all.IQs(data.all.idx==2,1),data.all.IQs(data.all.idx==2,2),blobBColor,...
        'MarkerSize',2,'DisplayName','Blob B')
    plot(data.all.C(:,1),data.all.C(:,2),'kx',...
        'MarkerSize',15,'LineWidth',3,'DisplayName','Fitted Centers')
    xlabel('I (V)')
    ylabel('Q (V)')
    title({'Clustered Quadrature Space Plot for ',...
        [num2str(length(data.gateI.Is)) ' Shots']})
    axis square
    grid on
    set(gca,'FontSize',14)
    plot([data.gateI.I data.gateX.I],[data.gateI.Q data.gateX.Q],'g+',...
        'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');

    legend('show')
    figure;
    histogram(data.all.IQsTransformed(:,1))
end


figure;
data.gateI.hist = histogram(data.gateI.IQsTransformed(:,2));
xValuesI = linspace(data.gateI.hist.BinLimits(1),data.gateI.hist.BinLimits(2),data.gateI.hist.NumBins);
fitI = fit(xValuesI', data.gateI.hist.Values', 'gauss1');             
hold on
plot(fitI, '-b');
hold on;
data.gateX.hist = histogram(data.gateX.IQsTransformed(:,2));
hold on;
xValuesX = linspace(data.gateX.hist.BinLimits(1),data.gateX.hist.BinLimits(2),data.gateX.hist.NumBins);
fitX = fit(xValuesX', data.gateX.hist.Values', 'gauss1');
plot(fitX, '-r');

twoSubGauss = @(x) fitX(x) - fitI(x);
xGuess = mean([fitX.b1 fitI.b1]);
xzero = fzero(twoSubGauss,xGuess);

%% Calculate the state preparation fidelities

% Assign the fitted functions to actual callable functions
funX = @(x) fitX(x);
funI = @(x) fitI(x);

% integrate over the entire span for both I and X
totalX = integral(funX,-inf,inf,'ArrayValued',true);
totalI = integral(funI,-inf,inf,'ArrayValued',true);

% integrate from I-X intersection point to INF for X gate
partialX = integral(funX,xzero,inf,'ArrayValued',true);

% integrate from -INF to I-X intersection point for I gate
partialI = integral(funI,-inf,xzero,'ArrayValued',true);

% Fidelity goes to 100% as intersection point moves out of actual gaussian
% regions
fidelX = partialX / totalX;
fidelI = partialI / totalI;

% display calculated fideilities on double gaussian plot
txtX = ['P_{|1>} = ',num2str(fidelX)];
txtI = ['P_{|0>} = ',num2str(fidelI)];
text(min(data.gateI.hist.BinLimits),fitX(xzero),txtX);
text(min(data.gateI.hist.BinLimits),fitX(xzero) + fitI(xzero)/10,txtI);

%% Calculate SNR from extracted gaussian parameters

% distance between the gaussian peaks
d = abs(fitX.b1 - fitI.b1);

% the factor of sqrt(2) comes in bc matlab defines its gaussian as 
% exp(- [(x - b)/c]^2)
sigmaX_I = fitI.c1 / sqrt(2);
sigmaX_X = fitX.c1 / sqrt(2);
SNR = d / sqrt(sigmaX_I^2 + sigmaX_X^2);
SNRdB = 10*log10(SNR);
txtSNR = ['SNR = ', num2str(SNRdB)];
text(min(data.gateI.hist.BinLimits),fitX(xzero) - fitX(xzero)/10,txtSNR);
end