function [data] = plotIQBlobsP1(CAL_I, CAL_X, DATA)
%% Select the calibration and data files
%calibration
if ~exist('CAL_I','var')
    [filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the file corresponding to the Identity calibration...',...
    'Select the file corresponding to the Pi calibration...'});
    
    if ~status
        return
    end
    
    calFileI = fullfile(pathnames{1}, filenames{1});
    calFileX = fullfile(pathnames{2}, filenames{2});
    calI = processMeasurementData(importMeasurementData(calFileI));
    calX = processMeasurementData(importMeasurementData(calFileX));
else
    calI = CAL_I;
    calX = CAL_X;
end

%data
if ~exist('DATA','var')
    data = loadMeasurementData;
else
    data = DATA;
end

% Call external function to find the threshold bisector and offset from the
% calibration files for |0> and |1>. rotInfo contains a shift and rotation
% matrix to threshold the data in the same fashion as the calibration files
[~, ~, ~, rotInfo] = fit2MeasIQBlobsToGaussHist(calI,calX);

%% Calculate P1 based on thresholding
% In the case of only one set of Is/Qs (no sweep parameter), only calculate
% P1 for the one set of reps
if length(data.Is(1,:)) == 1
    data.allIQs(:,:) = [data.Is(:) data.Qs(:)];
    data.allIQsTransformed(:,:) = shiftAndSpin(rotInfo, data.allIQs(:,:));
    if rotInfo.ZeroStateLoc > 0
        data.stateZero = length(data.allIQsTransformed(data.allIQsTransformed(:,1) > 0,1));
        data.stateOne = length(data.allIQsTransformed(data.allIQsTransformed(:,1) < 0,1));
    else
        data.stateZero = length(data.allIQsTransformed(data.allIQsTransformed(:,1) < 0,1));
        data.stateOne = length(data.allIQsTransformed(data.allIQsTransformed(:,1) > 0,1));
    end
    
        data.P1 = data.stateOne / (data.stateZero + data.stateOne);
% Calculate P1 as a function of a single sweep parameter.
else
    for n = 1:length(data.Is(:,1))
        data.allIQs(n,:,:) = [data.Is(n,:); data.Qs(n,:)]';
        data.allIQsTransformed(n,:,:) = shiftAndSpin(rotInfo, squeeze(data.allIQs(n,:,:)));
        if rotInfo.ZeroStateLoc > 0
            data.stateZero = length(data.allIQsTransformed(n,data.allIQsTransformed(n,:,1) > 0,1));
            data.stateOne = length(data.allIQsTransformed(n,data.allIQsTransformed(n,:,1) < 0,1));
        else
            data.stateZero = length(data.allIQsTransformed(n,data.allIQsTransformed(n,:,1) < 0,1));
            data.stateOne = length(data.allIQsTransformed(n,data.allIQsTransformed(n,:,1) > 0,1));
        end
    
        data.P1(n) = data.stateOne / (data.stateZero + data.stateOne);
    end
end

figure(); plot(data.(data.indep{1}),data.P1);
grid on
ylabel('P1')
set(gca,'FontSize',16)
xlabel([data.indep{1} ' (' data.units.(data.indep{1}) ')'])
end


function [ IQsTransformed ] = shiftAndSpin(rotInfo, IQs)
    IQsShifted = [(IQs(:,1)- rotInfo.ShiftI) (IQs(:,2) - rotInfo.ShiftQ)];
    IQsTransformed = IQsShifted * rotInfo.RMatrix;
end