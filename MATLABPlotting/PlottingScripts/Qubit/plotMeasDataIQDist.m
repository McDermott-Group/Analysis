function plotMeasDataIQDist
%plotMeasDataIQDist Plot IQ-distance between two data sets.

% Select files to compute the difference.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {'Select the first data file...',...
     'Select the second data file...'});
if ~status
    return
end

% Read the data files, convert the variable names, and specify the units.
file1 = fullfile(pathnames{1}, filenames{1});
file2 = fullfile(pathnames{2}, filenames{2});
data1 = processMeasurementData(importMeasurementData(file1));
data2 = processMeasurementData(importMeasurementData(file2));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, base_filename1] = fileparts(filenames{1});
[~, base_filename2] = fileparts(filenames{2});

for data_index = 1:length(data1.dep)
    I_name = data1.dep{data_index};

    if isempty(strfind(I_name, '_Std_Dev')) &&...
            ~isempty(strfind(I_name, 'I'))
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data1, Q_name)
            continue
        end
    else
        continue
    end

    if ~isfield(data2, I_name) && ~isfield(data2, Q_name)
        error('The selected files do not match.')
    else
        I1 = data1.(I_name);
        Q1 = data1.(Q_name);
        I2 = data2.(I_name);
        Q2 = data2.(Q_name);
    end
    I1_rels = data1.rels.(I_name);
    Q1_rels = data1.rels.(Q_name);
    I2_rels = data2.rels.(I_name);
    Q2_rels = data2.rels.(Q_name);
    
    if isempty(I1_rels) || isempty(Q1_rels) ||...
        isempty(I2_rels) || isempty(Q2_rels)
        continue
    end

    if length(I1_rels) ~= length(I2_rels) ||...
            length(Q1_rels) ~= length(I2_rels) ||...
            length(I2_rels) ~= length(Q1_rels)
        error('The selected files do not match.')
    end
    for k = 1:length(I1_rels)
        if ~strcmp(I1_rels{k}, I2_rels{k}) ||...
               ~strcmp(Q1_rels{k}, Q2_rels{k}) ||...
               ~strcmp(I1_rels{k}, Q1_rels{k}) ||...
               length(data1.(I1_rels{k})) ~= length(data2.(I2_rels{k})) ||...
               length(data1.(Q1_rels{k})) ~= length(data2.(Q2_rels{k})) ||...
               length(data1.(I1_rels{k})) ~= length(data2.(Q1_rels{k}))
            error('The selected files do not match.')
        end
    end
    
    for k = 1:length(I1_rels)
        indep1 = data1.(I1_rels{k});
        indep2 = data2.(I2_rels{k});
        if any(indep1 ~= indep2)
            if all(flip(indep2) == indep1)
                disp(I2_rels{k})
                I2 = flip(I2, k);
                Q2 = flip(Q2, k);
            else
                error('The selected files do not match.')
            end
        end
    end
    
    try
        distance = sqrt((I2 - I1).^2 + (Q2 - Q1).^2);
    catch
        error('The selected files do not match.')
    end

    data = data1;
    dist_name = strrep(I_name, 'I', 'Quadrature_Distance');
    data.(dist_name) = distance;
    data.units.(dist_name) = data1.units.(I_name);
    data.rels.(dist_name) = I1_rels;
    data.dep{length(data.dep) + 1} = dist_name;
    data.plotting.(dist_name).full_name = strrep(dist_name, '_', ' ');
    data.plotting.(dist_name).plot_title =...
        {[strrep(dist_name, '_', ' '), getUnits(data1, I_name),...
        ' between Two Datasets:'],...
         [strrep(filenames{1}, '_', '\_'), ' [', data1.Timestamp, ']'],...
         [strrep(filenames{2}, '_', '\_'), ' [', data2.Timestamp, ']']};
    data.plotting.(dist_name).plot_filename =...
        fullfile(plts_path, [base_filename1, '-', base_filename2,...
        '_', dist_name]);
     
    plotDataVar(data, dist_name);
end