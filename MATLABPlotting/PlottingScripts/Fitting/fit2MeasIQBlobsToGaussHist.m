function [data, maxFidelity, singleShotFidelity, probOne, rotInfo] = ...
        fit2MeasIQBlobsToGaussHist(GATE_I, GATE_X, PLOT_HIST)
% This function is meant to take data from an identitiy gate and a Pi gate
% and find the centroids of the IQ space blobs, center and rotate them
% along the I axis and take a line cut through them. Additionally, it will
% integrate the data and plot the results ERF's of the gaussian
% distribution.

% Number of bins in a histogram.
Nbins = 250;

% Check to see if gate data was fed into the function or needs to be.
if ~exist('GATE_I','var')
    % Select files to compute the difference.
    [filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the file corresponding to the Identity gate...',...
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

% Initialize the data stack
data.gateI.Is = gateI.Is();
data.gateI.Qs = gateI.Qs();
data.gateX.Is = gateX.Is();
data.gateX.Qs = gateX.Qs();
 
% Combine gateI and gateX I quadrature data.
data.allIs = [data.gateI.Is; data.gateX.Is];

% Combine gateI and gateX Q quadrature data.
data.allQs = [data.gateI.Qs; data.gateX.Qs];

% Combine all I-Q data into one array where the I gate data is the 1st half 
% (1:end/2) and X Gate data is the 2nd half of the array (end/2 : end).
data.allIQs = [data.allIs data.allQs];

% Call centerAndRoate function that finds the IQ blob centers of
% the combined I and X data sets
[data.allIQsRot, rotInfo] = centerAndRotate(data.allIQs);

% Histogram the data and create requisite arrays for fitting
[Icounts, Iedges] = histcounts(data.allIQsRot(1:end/2,1), Nbins);
[Xcounts, Xedges] = histcounts(data.allIQsRot(end/2:end,1), Nbins);

xValuesI = zeros(length(Iedges)-1,1);
xValuesX = zeros(length(Xedges)-1,1);
for i=1:length(Iedges)-1
   xValuesI(i) = mean([Iedges(i) Iedges(i+1)]);
   xValuesX(i) = mean([Xedges(i) Xedges(i+1)]);
end

set1 =(data.allIQsRot(1:end/2,1));
set2 = (data.allIQsRot(end/2+1:end,1));
fulldata = [set1;set2];

[fullcounts, fulledges] = histcounts(fulldata, Nbins);
fullValues = zeros(length(fulledges)-1,1);

for i=1:length(fulledges)-1
   fullValues(i) = mean([fulledges(i) fulledges(i+1)]);
end

fitfull = fit(fullValues, fullcounts',  'gauss2');

a1=fitfull.a1;
b1=fitfull.b1;
c1=fitfull.c1;

a2=fitfull.a2;
b2=fitfull.b2;
c2=fitfull.c2;

% Fit the data and determine the location of the "1" state

fitI = fit(xValuesI,Icounts','gauss1');
fitX = fit(xValuesX,Xcounts','gauss1');


fittwogaussianX = a1*exp(-((xValuesX-b1)/c1).^2);
fittwogaussianI = a2*exp(-((xValuesI-b2)/c2).^2);


if fitfull.b1 > 0
    rotInfo.ZeroStateLoc = 1;
else
    rotInfo.ZeroStateLoc = -1;
end

% Integrate over the histograms to learn the cumulative probatility of each
% given outcome. 

intX = cumtrapz(fittwogaussianX);
intI = cumtrapz(fittwogaussianI);

funX = @(x) fitX(x);

probOne = integral(funX, -inf, xValuesX(end/2), 'ArrayValued', true) /...
    integral(funX, -inf, inf, 'ArrayValued', true);

% Use the fitted values of the gaussian histograms to find the maximum 
% measurement fidelity (separation fidelity)
fidelity =  @(x) 0.5 * abs((erf((x-fitfull.b1)/(fitfull.c1))) - (erf((x-fitfull.b2)/(fitfull.c2))));
allX = [xValuesX xValuesI];
fids = zeros(size(allX));
for n = 1:length(allX)
    fids(n) = fidelity(allX(n));
end
maxFidelity = max(fids(:));

% Interpolate one of the cumtrapz so the indicies match up with the other
% for a fair comparison. This will ensure that doing simple subtraction via
% indicies actually compares like X values

intX_interp = interp1(xValuesX,intX,xValuesI);


[singleShotFidelity, SSFMax_idx] = max(abs(intX_interp/max(intX) - ...
                                           intI/max(intI)));


if plotHist
    createFigure([.9, .1, .88, .8]);
    histogram(data.allIQsRot(1:end/2,1), Nbins); 
    hold on
    histogram(data.allIQsRot(end/2:end,1), Nbins);
    hold on
    histogram(fulldata,Nbins)

    

    h = plot(xValuesI,fittwogaussianI);
    set(h, 'LineWidth',2);   
    h = plot(xValuesX,fittwogaussianX);
    set(h, 'LineWidth',2);     
    


    
    ylabel('Counts', 'FontSize', 14)
    % Note: yyaxis not available before 2016a
    try
        yyaxis right
        h = plot(xValuesI, intI/max(intI), '--b');
        ylabel('Occupation Probability', 'FontSize', 14)
        ax = gca;
        ax.YColor = [0,0,0];
        set(h, 'LineWidth',2);
        h = plot(xValuesX, intX/max(intX), '--r');
        set(h, 'LineWidth',2);
        line([xValuesI(SSFMax_idx), xValuesI(SSFMax_idx)],...
             [intI(SSFMax_idx)/max(intI), intX_interp(SSFMax_idx)/max(intX)],...
             'DisplayName', 'Max Fidelity',...
             'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)
    catch
    end
    
    grid on
    xunits = getUnits(gateI, 'Is');
    xlabel(['Generalized Quadrature Coordinate', xunits], 'FontSize', 14)
    [~, filenameI, extI] = fileparts(gateI.Filename);
    [~, filenameX, extX] = fileparts(gateX.Filename);
    title({['Separation Fidelity = ',...
                num2str(100 * maxFidelity, '%.3f'), '%'],...
           ['Single Shot Fidelity = ',...
                num2str(100 * singleShotFidelity, '%.3f'), '%'],...
           ['Dataset A: ', strrep(filenameI, '_', '\_'), extI,...
            ' [', gateI.Timestamp, ']'],...
           ['Dataset B: ', strrep(filenameX, '_', '\_'), extX,...
            ' [', gateX.Timestamp, ']']}, 'FontSize', 14)
     legend('dataset A', 'dataset B',...
         'gaussian fit to A', 'gaussian fit to B',...
         'integrated dataset A', 'integrated dataset B',...
         'Location', 'NorthWest')
end
end

function [IQshiftRot, rotInfo] = centerAndRotate(allIQs)
%centerAndRotate(allIQs) is designed to take two IQ blobs, find the center
%  between them, shift that center to zero, and then rotate the two blobs
%  such that the Q=0 line bisects the two blobs. The function then returns
%  this information and the rotation matrix such that other functions may
%  perform the exact same rotation and shift in IQ space.

    opts = statset('UseParallel', 0,...
                   'MaxIter', 10000);

    [~, C] = kmeans(allIQs, 2, 'Replicates', 1, 'Options', opts);

    % Find the center between the 2 blobs
    rotCenterI = mean(C(:,1));
    rotCenterQ = mean(C(:,2));

    IQshift = [allIQs(:,1)- rotCenterI, allIQs(:,2) - rotCenterQ];
    centShift = [C(:,1) - rotCenterI, C(:,2) - rotCenterQ];
    angle = atan(centShift(1,2) / centShift(1,1));

    rotMatrix = [cos(angle) -sin(angle);...
                 sin(angle)  cos(angle)]; 

    IQshiftRot = IQshift * rotMatrix;

    rotInfo.RMatrix = rotMatrix;
    rotInfo.ShiftI = rotCenterI;
    rotInfo.ShiftQ = rotCenterQ;
end

