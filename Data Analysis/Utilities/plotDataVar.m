function plotDataVar(data, dependent_variable, type)
%plotDataVar    Create plots for a data variable.
%
%   plotDataVar(DATA, DEPENDENT_VARIABLE) creates plots for
%   DEPENDENT_VARIABLE in structure DATA. Specify plot TYPE if only
%   a specific plot is desired. Supported plot types: 'simple' (1D),
%   'errorbar' (1D), 'smooth' (2D), 'pixelated' (2D), and 'polar' (2D).

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
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([full_dep_name, yunits], 'FontSize', 14);
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
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14);
        ylabel([full_dep_name, yunits],...
            'FontSize', 14);
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
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits, ':'],...
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
            title({[full_dep_name, zunits, ':'],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']'], indep_vars}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_smooth']);
    end
    if strcmp(type, '') || strcmp(type, 'pixelated')
        % Plot the data as a pixelated image.
        createFigure('right');
        plotPixelated(indep_vals1, indep_vals2, dep_vals');
        xlabel([strrep(indep_name1, '_', ' '), xunits], 'FontSize', 14);
        ylabel([strrep(indep_name2, '_', ' '), yunits], 'FontSize', 14);
        if exist('plot_title', 'var')
            title(plot_title, 'FontSize', 10)
        else
            title({[full_dep_name, zunits, ':'],...
                   [strrep(filename, '_', '\_'), ext,...
                   ' [', data.Timestamp, ']']}, 'FontSize', 10)
        end
        savePlot([plot_filename, extra_filename, '_pixelated']);
    end
end
if length(dep_rels) > 2 && strcmp(type, '')
    disp(['Data variable ''', strrep(dependent_variable, '_', ' '),...
          ''' depends on more than two sweep variables. ',...
          'This data will not be plotted.'])
end