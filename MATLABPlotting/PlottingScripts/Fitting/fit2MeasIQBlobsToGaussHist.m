function [data, maxFidelity, probOne, rotInfo] = fit2MeasIQBlobsToGaussHist(GATE_I, GATE_X, PLOT_HIST)
%% This function is meant to take data from an identitiy gate and a Pi gate
% and find the centroids of the IQ space blobs, center and rotate them
% along the I axis and take a line cut through them. Additionally, it will
% integrate the data and plot the results ERF's of the gaussian
% distribution

%% Check to see if gate data was fed into the function or needs to be
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

if ~exist('PLOT_HIST','var')
    plotHist = 1;
else
    plotHist = PLOT_HIST;
end

%% Initialize the data stack
data.gateI.Is = gateI.Is();
data.gateI.Qs = gateI.Qs();
data.gateX.Is = gateX.Is();
data.gateX.Qs = gateX.Qs();
 
% Combine gateI and gateX I quad data
data.allIs = [data.gateI.Is; data.gateX.Is];

% Combine gateI and gateX Q quad data
data.allQs = [data.gateI.Qs; data.gateX.Qs];

% combine all I-Q data into one array where the I gate data is 1st half (1:end/2) and X Gate data is the 2nd half
% of the array (end/2 : end)
data.allIQs = [data.allIs data.allQs];

%% Call centerAndRoate function that finds the IQ blob centers of the combined I and X data sets
[data.allIQsRot, rotInfo] = centerAndRotate(data.allIQs);


if plotHist
    
    figure;
    histI = histogram(data.allIQsRot(1:end/2,1),250); hold on;
    histX = histogram(data.allIQsRot(end/2:end,1),250);

    xValuesI = linspace(histI.BinLimits(1),histI.BinLimits(2),histI.NumBins);
    fitI = fit(xValuesI', histI.Values', 'gauss1');
    if fitI.b1 > 0
        rotInfo.ZeroStateLoc = 1;
    else
        rotInfo.ZeroStateLoc = -1;
    end

    xValuesX = linspace(histX.BinLimits(1),histX.BinLimits(2),histX.NumBins);
    %options = fitoptions('gauss2', 'StartPoint', [100 -fitI.b1 1 100 fitI.b1 1]); 
    fitX = fit(xValuesX', histX.Values', 'gauss2');

    h = plot(fitI, '-b');
    set(h, 'LineWidth',2);
    h = plot(fitX, '-r');
    set(h, 'LineWidth',2);

    funX = @(x) fitX(x);
    
    intX = cumtrapz(xValuesX,histX.Values);
    probOne = integral(funX,-inf,xValuesX(end/2),'ArrayValued',true)/integral(funX,-inf,inf,'ArrayValued',true);
    intI = cumtrapz(xValuesI,histI.Values);
    
    yyaxis right;
    h = plot(xValuesI,intI/max(intI), '-b');
    ax = gca;
    ax.YColor = [0,0,0];
    set(h, 'LineWidth',2);
    h = plot(xValuesX,intX/max(intX), '-r');
    set(h, 'LineWidth',2);
    fidelity = intI/max(intI) - intX/max(intX);
    fidelity = abs(fidelity);
    maxFidelity = max(fidelity);
    
else
    
    [Icounts, Iedges] = histcounts(data.allIQsRot(1:end/2,1),250);
    % [Xcounts, Xedges] = histcounts(data.allIQsRot(end/2:end,1),250);
    xValuesI = zeros(length(Iedges)-1,1);
    for i=1:length(Iedges)-1
       xValuesI(i) = mean([Iedges(i) Iedges(i+1)]);
       % xValuesX = mean([Xedges(i) Xedges(i+1)]);
    end
    fitI = fit(xValuesI, Icounts', 'gauss1');
    if fitI.b1 > 0
        rotInfo.ZeroStateLoc = 1;
    else
        rotInfo.ZeroStateLoc = -1;
    end
    maxFidelity = Inf;
    probOne = Inf;
end

end

function [IQshiftRot, rotInfo] = centerAndRotate(allIQs)
%centerAndRotate(allIQs) is designed to take two IQ blobs, find the center
%  between them, shift that center to zero, and then rotate the two blobs
%  such that the Q=0 line bisects the two blobs. The function then returns
%  this information and the rotation matrix such that other functions may
%  perform the exact same rotation and shift in IQ space.

opts = statset('Display','final',...
               'UseParallel', 0,...
               'MaxIter', 10000);

[~, C] = kmeans(allIQs, 2, 'Replicates', 1,...
                                  'Options', opts);

% Find the center between the 2 blobs
rotCenterI = mean(C(:,1));
rotCenterQ = mean(C(:,2));

IQshift = [allIQs(:,1)- rotCenterI allIQs(:,2) - rotCenterQ];
centShift = [C(:,1) - rotCenterI C(:,2) - rotCenterQ];
angle = atan(centShift(1,2)/centShift(1,1));

rotMatrix = [cos(angle) -sin(angle);...
             sin(angle)  cos(angle)]; 
IQshiftRot = IQshift * rotMatrix;

rotInfo.RMatrix = rotMatrix;
rotInfo.ShiftI = rotCenterI;
rotInfo.ShiftQ = rotCenterQ;

end

