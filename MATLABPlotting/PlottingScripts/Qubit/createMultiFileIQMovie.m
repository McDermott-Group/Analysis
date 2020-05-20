function createIQMovie
%createIQMovie   Animate the I-Q space blob evolution.


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

meanI = [];
meanQ = [];
center0I = [];
center0Q = [];
center1I = [];
center1Q = [];
j = 0;

for i = 1:length(filenames)
    try
        file = fullfile(pathnames{i}, filenames{i});
        data = loadMeasurementData(file);
    catch
        warning(['The file ' filenames{i} ' contains unmatched data.'])
        continue
    end
    
    if ~exist('I_name', 'var')
        I_name = selectDepDataVars(data);
        I_name = I_name{1};
        Q_name = selectDepDataVars(data);
        Q_name = Q_name{1};
    end
    
    if strcmp(data.rels.(I_name){1}, 'Repetition_Index') ||...
            strcmp(data.rels.(I_name){1}, 'Long_Repetition_Index')
        I = data.(I_name);
        Q = data.(Q_name);
        if length( data.rels.(I_name) ) > 1
            indep = data.rels.(I_name){2};
        end
    elseif strcmp(data.rels.(I_name){2}, 'Repetition_Index') ||...
            strcmp(data.rels.(I_name){2}, 'Long_Repetition_Index')
        I = data.(I_name)';
        Q = data.(Q_name)';
        if length( data.rels.(I_name) ) > 1
            indep = data.rels.(I_name){1};
        end
    else
        continue
    end
    
    if length( data.rels.(I_name) ) > 1
        indep_vals = data.(indep);
        iunits = getUnits(data, indep);
    else
        indep_vals = [1];
        iunits = '';
        indep = '';
    end
    meanI = [meanI mean(I)];
    meanQ = [meanQ mean(Q)];
    if isfield(data, 'Center_of_0_State_I')
        center0I = [center0I data.Center_of_0_State_I];
        center0Q = [center0Q data.Center_of_0_State_Q];
    else
        center0I = [center0I nan];
        center0Q = [center0Q nan];
    end
    if isfield(data, 'Center_of_1_State_I')
        center1I = [center1I data.Center_of_1_State_I];
        center1Q = [center1Q data.Center_of_1_State_Q];
    else
        center1I = [center1I nan];
        center1Q = [center1Q nan];
    end
    
    if i == 1
        % Create folder Plots if necessary.
        mov_path = makeDirPlots(pathnames{i}, 'Movies');
        % Create a new video file.
        vidObj = VideoWriter(fullfile(mov_path, [filenames{i}, '_', I_name, '-',...
                Q_name, '.avi']));
        vidObj.FrameRate = 15; %ceil(length(indep_vals) / 15);
        open(vidObj);

        createFigure;
        xunits = getUnits(data, I_name);
        yunits = getUnits(data, Q_name);

        xl = [min(I(:)) max(I(:))];
        yl = [min(Q(:)) max(Q(:))];
        plot(xl, yl)
        xlim(xl)
        ylim(yl)
        axis equal
        axs = axis;
    end
    
    for k = 1:length(indep_vals)
        j = j + 1;
        if size(I,1) > 1
            plot(I(:, k), Q(:, k), '.', 'MarkerSize', 10)
        else
            plot(I, Q, '.', 'MarkerSize', 10)
        end
        hold on
        plot(meanI(1:j), meanQ(1:j), 'r-', 'LineWidth', 2)
        plot(meanI(j), meanQ(j), 'k+', 'MarkerSize', 10)
        plot(center0I(1:j), center0Q(1:j), 'g-', 'LineWidth', 2)
        plot(center0I(j), center0Q(j), 'k+', 'MarkerSize', 10)
        plot(center1I(1:j), center1Q(1:j), 'c-', 'LineWidth', 2)
        plot(center1I(j), center1Q(j), 'k+', 'MarkerSize', 10)
        hold off
        axis(axs);
        grid on
        title({[filenames{i}, ' [', data.Timestamp, ']'],...
            [strrep(indep, '_', ' '), ' = ', num2str(indep_vals(k)), iunits]},...
            'Interpreter', 'none', 'FontSize', 10)
        xlabel([strrep(I_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(Q_name, '_', ' '), yunits], 'FontSize', 14)
        % Write each frame to the file.
        writeVideo(vidObj, getframe(gcf));
    end
    
    if i == length(filenames)
        % Close the file.
        close(vidObj);
    end
end
