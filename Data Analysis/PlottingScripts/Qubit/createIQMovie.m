function createIQMovie
%createIQMovie   Animate the I-Q space blob evolution.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    
    if ~isempty(strfind(I_name, '_Std_Dev')) || isempty(strfind(I_name, 'I'))
        continue
    else
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data, Q_name) ||...
                length(data.rels.(I_name)) ~= 2 ||...
                length(data.rels.(Q_name)) ~= 2 ||...
                ~strcmp(data.rels.(I_name){1}, data.rels.(Q_name){1}) ||...
                ~strcmp(data.rels.(I_name){2}, data.rels.(Q_name){2}) 
            continue
        end
    end
    
    if strcmp(data.rels.(I_name){1}, 'Repetition_Index') ||...
            strcmp(data.rels.(I_name){1}, 'Long_Repetition_Index')
        I = data.(I_name);
        Q = data.(Q_name);
        indep = data.rels.(I_name){2};
    elseif strcmp(data.rels.(I_name){2}, 'Repetition_Index') ||...
            strcmp(data.rels.(I_name){2}, 'Long_Repetition_Index')
        I = data.(I_name)';
        Q = data.(Q_name)';
        indep = data.rels.(I_name){1};
    else
        continue
    end
    
    indep_vals = data.(indep);
    meanI = mean(I);
    meanQ = mean(Q);
    
    % Create folder Plots if necessary.
    mov_path = makeDirPlots(pathname, 'Movies');
    % Create a new video file.
    vidObj = VideoWriter(fullfile(mov_path, [filename, '_', I_name, '-',...
            Q_name, '.avi']));
    vidObj.FrameRate = ceil(length(indep_vals) / 15);
    open(vidObj);

    createFigure;
    xunits = getUnits(data, I_name);
    yunits = getUnits(data, Q_name);
    iunits = getUnits(data, indep);
    xlabel([strrep(I_name, '_', ' '), xunits], 'FontSize', 14);
    ylabel([strrep(Q_name, '_', ' '), yunits], 'FontSize', 14);

    axis([min(I(:)) max(I(:)) min(Q(:)) max(Q(:))])
    axis equal
    axs = axis;
    for k = 1:length(indep_vals)
        plot(I(:, k), Q(:, k), '.', 'MarkerSize', 10)
        hold on
            plot(meanI(1:k), meanQ(1:k), 'r-', 'LineWidth', 2)
            plot(meanI(k), meanQ(k), 'k+', 'MarkerSize', 10)
        hold off
        axis(axs)
        grid on
        title({[filename, ext, ' [', data.Timestamp, ']'],...
            [strrep(indep, '_', ' '), ' = ', num2str(indep_vals(k)), iunits]},...
            'Interpreter', 'none', 'FontSize', 10)
        % Write each frame to the file.
        writeVideo(vidObj, getframe(gcf));
    end
    % Close the file.
    close(vidObj);
end
end