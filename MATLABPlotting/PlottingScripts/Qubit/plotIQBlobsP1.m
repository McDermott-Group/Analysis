function plotIQBlobsP1(CAL_I, CAL_X, DATA, MAKE_PLOTS)

% Select the calibration files.
if ~exist('CAL_I', 'var')
    [filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the file corresponding to the Identity calibration...',...
     'Select the file corresponding to the Pi calibration...'});
    
    if ~status
        return
    end
    
    calFileI = fullfile(pathnames{1}, filenames{1});
    calFileX = fullfile(pathnames{2}, filenames{2});
    calI = loadMeasurementData(calFileI);
    calX = loadMeasurementData(calFileX);

    
else
    calI = CAL_I;
    calX = CAL_X;
end

% Select the data file.
if ~exist('DATA', 'var')
    data = loadMeasurementData;
else
    data = DATA;
end

if ~exist('MAKE_PLOTS', 'var')
    makePlots = 1;
else
    makePlots = MAKE_PLOTS;
end

% Call external function to find the threshold bisector and offset from the
% calibration files for |0> and |1>. rotInfo contains a shift and rotation
% matrix to threshold the data in the same fashion as the calibration files



%This patch allows us to plot .mat files, .mat files have [1x4000]
%matrices, instead of the 4000 x 1 matrices that txt and hdf5 data sets
%provide, this patch simply transposes the data that needs to be altered.

if size(calX.Qs,1) == 1

    calX.Qs = calX.Qs';
    calX.Repetition_Index = calX.Repetition_Index';
    calX.Amplitudes = calX.Amplitudes';
    calX.Is = calX.Is';
    calX.Phases = calX.Phases';

    calI.Qs = calI.Qs';
    calI.Amplitudes = calI.Amplitudes';
    calI.Is = calI.Is';
    calI.Phases = calI.Phases';

end


[~, ~, ~, ~, rotInfo] = fit2MeasIQBlobsToGaussHist(calI, calX, makePlots);

if isempty(fields(data))
    return
end
% Calculate P1 based on thresholding.
% In the case of only one set of Is/Qs (no sweep parameter), only calculate
% P1 for the one set of reps.
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
else % Calculate P1 as a function of a single sweep parameter.
    for n = 1:length(data.Is(:,1))
        data.allIQs(n,:,:) = [data.Is(n,:); data.Qs(n,:)]';
        data.allIQsTransformed(n,:,:) = shiftAndSpin(rotInfo, squeeze(data.allIQs(n,:,:)));
        if rotInfo.ZeroStateLoc > 0
            data.stateZero = length(data.allIQsTransformed(n, data.allIQsTransformed(n,:,1) > 0, 1));
            data.stateOne = length(data.allIQsTransformed(n, data.allIQsTransformed(n,:,1) < 0, 1));
        else
            data.stateZero = length(data.allIQsTransformed(n, data.allIQsTransformed(n,:,1) < 0, 1));
            data.stateOne = length(data.allIQsTransformed(n, data.allIQsTransformed(n,:,1) > 0, 1));
        end
    
        data.P1(n) = data.stateOne / (data.stateZero + data.stateOne);
    end
end

if makePlots
    occup_prob = 'Excited_State_Occupation_Probability';
    data.(occup_prob) = data.P1;

    data.units.(occup_prob) = '';
    data.rels.(occup_prob) = data.rels.I;
    data.dep{length(data.dep) + 1} = occup_prob;
    data.plotting.(occup_prob).full_name =...
        'Excited State Occupation Probability';
    [pathname, filename, ext] = fileparts(data.Filename);
    plts_path = makeDirPlots(pathname);
    [~, filenameI, extI] = fileparts(calI.Filename);
    [~, filenameX, extX] = fileparts(calX.Filename);
    data.plotting.(occup_prob).plot_title =...
        {['Data: ', strrep(filename, '_', '\_'), ext,...
        ' [', calI.Timestamp, ']'],...
        ['Calibration I: ', strrep(filenameI, '_', '\_'), extI,...
        ' [', calI.Timestamp, ']'],...
        ['Calibration X: ', strrep(filenameX, '_', '\_'), extX,...
        ' [', calX.Timestamp, ']']};
    data.plotting.(occup_prob).plot_filename =...
            fullfile(plts_path, [filename, '_occup_prob']);
    
    plotDataVar(data, occup_prob)
    ylim([0, 1])
end

end


function [ IQsTransformed ] = shiftAndSpin(rotInfo, IQs)
    IQsShifted = [(IQs(:,1)- rotInfo.ShiftI) (IQs(:,2) - rotInfo.ShiftQ)];
    IQsTransformed = IQsShifted * rotInfo.RMatrix;
end