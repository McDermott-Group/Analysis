function plotDataVar(data, dependent_variable, type)
%plotDataVar    Create plots for a data variable.
%
%   plotDataVar(DATA, DEPENDENT_VARIABLE) creates plots for
%   DEPENDENT_VARIABLE in structure DATA. Specify plot TYPE if only
%   a specific plot is desired. Supported plot types: 'simple' (1D),
%   'errorbar' (1D), 'smooth' (2D), 'pixelated' (2D), 'polar' (2D), and
%   'scatter' (3D scatter plot of 2D).

if isempty(fields(data))
    return
end

if ~exist('type', 'var')
    type = '';
else
    type = lower(type);
end
if ~isfield(data, dependent_variable)
    if isfield(data, strrep(dependent_variable, '_', ' '))
        dependent_variable = strrep(dependent_variable, '_', ' ');
    elseif isfield(data, strrep(dependent_variable, ' ', '_'))
        dependent_variable = strrep(dependent_variable, ' ', '_');
    else
        error(['Variable ''', dependent_variable,...
            ''' is not found in the data structure.'])
    end
end

[pathname, filename, ext] = fileparts(data.Filename);
plts_path = makeDirPlots(pathname);

full_dep_name = strrep(dependent_variable, '_', ' ');
plot_filename = fullfile(plts_path, [filename, '_', dependent_variable]);
extra_filename = '';
if isfield(data, 'plotting') && isfield(data.plotting, dependent_variable)
    if isfield(data.plotting.(dependent_variable), 'full_name')
        full_dep_name = data.plotting.(dependent_variable).full_name;
    end
    if isfield(data.plotting.(dependent_variable), 'plot_filename')
        plot_filename = data.plotting.(dependent_variable).plot_filename;
    end
    if isfield(data.plotting.(dependent_variable), 'extra_filename')
        extra_filename = data.plotting.(dependent_variable).extra_filename;
    end
    if isfield(data.plotting.(dependent_variable), 'plot_title')
        plot_title = data.plotting.(dependent_variable).plot_title;
    end
end

dep_vals = data.(dependent_variable);
dep_rels = data.rels.(dependent_variable);
if isempty(dep_rels)
    disp(['Independent (sweep) variables for data variable ''',...
          strrep(dependent_variable, '_', ' '), ''' are not specified.'])
end

% Plot 1D data.
if length(dep_rels) == 1
    indep_name = dep_rels{1};
    indep_vals = data.(indep_name);

    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, dependent_variable);

    % Plot an errobar graph.
    if (strcmp(type, '') || strcmp(type, 'errorbar')) &&...
            isfield(data, 'error') &&...
            isfield(data.error, dependent_variable)
        createFigure('right');
        err = data.error.(dependent_variable);
        if (size(err, 1) == 1 || size(err, 2) == 1) &&...
                length(err) == length(indep_vals)
            plotErrorbar(indep_vals, dep_vals, err)
        elseif size(err, 1) == 2 && size(err, 2) == length(indep_vals)
            plotAssymErrorbar(indep_vals, dep_vals, err(1, :), err(2, :))
        elseif size(err, 2) == 2 && size(err, 1) == length(indep_vals)
             plotAssymErrorbar(indep_vals, dep_vals, err(:, 1), err(:, 2))
        end
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([full_dep_name, yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title({strrep(plot_title{1}, yunits, ''),...
                plot_title{2:end}}, 'FontSize', 10)
        else
            title({[strrep(filename, '_', '\_'), ext,...
                ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_errorbar']);
    end

    % Plot a simple 1D graph.
    if strcmp(type, '') || strcmp(type, 'simple')
        createFigure;
        plotSimple(indep_vals, dep_vals)
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([full_dep_name, yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title({strrep(plot_title{1}, yunits, ''),...
                plot_title{2:end}}, 'FontSize', 10)
        else
            title({[strrep(filename, '_', '\_'), ext,...
                ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_simple']);
    end
end

% Plot 2D data.
if length(dep_rels) == 2
    indep_name1 = dep_rels{1};
    indep_name2 = dep_rels{2};
    indep_vals1 = data.(indep_name1);
    indep_vals2 = data.(indep_name2);

    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    zunits = getUnits(data, dependent_variable);

    if (strcmp(type, '') || strcmp(type, 'smooth')) &&...
            isempty(strfind(indep_name1, 'Phase')) &&...
            isempty(strfind(indep_name2, 'Phase'))
        % Plot the data as a smooth surface.
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_smooth']);
    elseif strcmp(type, '') || strcmp(type, 'polar')
        % Create a polar (smooth) plot.
        if ~isempty(strfind(indep_name1, 'Phase'))
            phase = indep_vals1;
            radius = indep_vals2;
            vals = dep_vals';
            indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                yunits, '; Phase: ', strrep(indep_name1, '_', ' '), xunits];
        else
            phase = indep_vals2;
            radius = indep_vals1;
            vals = dep_vals;
            indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                xunits, '; Phase: ', strrep(indep_name2, '_', ' '), yunits];
        end
        createFigure;
        plotPolar(radius, phase, vals);
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']'], indep_vars}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_smooth']);
    end
    if strcmp(type, '') || strcmp(type, 'pixelated')
        % Plot the data as a pixelated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_pixelated']);
    end
    if strcmp(type, 'scatter')
        createFigure('left');
        plot3DSimple(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_3d_scatter']);
    end
end



if length(dep_rels) == 3
%     freqIndexArray = size(dep_vals(1,1,:));
%     freqIndexArray = freqIndexArray(3);
%     for ii = 1:freqIndexArray
%         ff = BiasedLorentzian(data.RF_Frequency, reshape(dep_vals(ii,ii,:),[freqIndexArray,1]));
%     end
        
    dep_vals = min(dep_vals,[],3);
    [dep_vals1, index_min_s21] = min(dep_vals,[],3);
    indep_name1 = dep_rels{1};
    indep_name2 = dep_rels{2};
    indep_vals1 = data.(indep_name1);
    indep_vals2 = data.(indep_name2);

    xunits = getUnits(data, indep_name1);
    yunits = getUnits(data, indep_name2);
    zunits = getUnits(data, dependent_variable);

    if (strcmp(type, '') || strcmp(type, 'smooth')) &&...
            isempty(strfind(indep_name1, 'Phase')) &&...
            isempty(strfind(indep_name2, 'Phase'))
        % Plot the data as a smooth surface.
        createFigure;
        plotSmooth(indep_vals1, indep_vals2, dep_vals);
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_smooth']);
    elseif strcmp(type, '') || strcmp(type, 'polar')
        % Create a polar (smooth) plot.
        if ~isempty(strfind(indep_name1, 'Phase'))
            phase = indep_vals1;
            radius = indep_vals2;
            vals = dep_vals';
            indep_vars = ['Radius: ', strrep(indep_name2, '_', ' '),...
                yunits, '; Phase: ', strrep(indep_name1, '_', ' '), xunits];
        else
            phase = indep_vals2;
            radius = indep_vals1;
            vals = dep_vals;
            indep_vars = ['Radius: ', strrep(indep_name1, '_', ' '),...
                xunits, '; Phase: ', strrep(indep_name2, '_', ' '), yunits];
        end
        createFigure;
        plotPolar(radius, phase, vals);
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']'], indep_vars}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_smooth']);
    end
    if strcmp(type, '') || strcmp(type, 'pixelated')
        % Plot the data as a pixelated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_pixelated']);
    end
    if strcmp(type, 'scatter')
        createFigure('left');
        plot3DSimple(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14)
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_3d_scatter']);
    end
    
%     data_variable
%     processed_data_var = ['CutSubtracted_', data_variable];
%     data.(processed_data_var) = dep_vals;
%     data.units.(processed_data_var) = data.units.(data_variable);
%     data.rels.(processed_data_var) = data.rels.(data_variable);
%     data.dep{length(data.dep) + 1} = processed_data_var;
%     data.plotting.(processed_data_var).full_name =...
%         ['Cut-Subtracted ', strrep(data_variable, '_', ' ')];
%     data.plotting.(processed_data_var).extra_filename =...
%         ['_', normalization_direction];
% 
%     plotDataVar(data, processed_data_var);
%     
%     saveMeasData(data, [filename, '_', data_variable, '_cut_subtr'])
end


if length(dep_rels) > 3 && strcmp(type, '')
    disp(['Data variable ''', strrep(dependent_variable, '_', ' '),...
          ''' depends on more than three sweep variables. ',...
          'This data will not be plotted.'])
end
end 






function f = BiasedLorentzian(x, y)

median_y = median(y);
max_y = max(y);
min_y = min(y);
if max_y - median_y >= median_y - min_y
    background = min_y;
    idx = find(y == max_y, 1, 'first');
else
    background = max_y;
    idx = find(y == min_y, 1, 'first');
end

start = x(idx);
half_y = (median_y + y(idx)) / 2;
idx1 = find(y < half_y, 1, 'first');
idx2 = find(y < half_y, 1, 'last');
if ~isempty(idx1) && ~isempty(idx2) && x(idx1) ~= x(idx2)
    width = .05 * abs(x(idx2) - x(idx1));
else
    width = 0.01;
end

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.TolX = 1e-20;
opts.TolFun = 1e-20;
opts.StartPoint = [.5 * (max_y - min_y) * width^2 / 4,...
    start, width, 0, background];

f = fit(x(:), y(:), '(a / ((x - b)^2 + (c / 2)^2)) + d * (x - b) + e', opts);

format long
f.b

end