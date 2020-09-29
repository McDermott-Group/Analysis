function fitMeasData2QPT1Decay(T1_expt, data_variable, data)
%fitMeasData2QPT1Decayu(T1_EXPT, DATA_VARIABLE, DATA)  Fit T1 data to a
%multi-part exponential function with ability to feed in a known T1

%   fitMeasData2Exponent(T1_expt, DATA_VARIABLE, DATA) fits data for
%   DATA_VARIABLE to an exponent with a Poisson distribution of
%   quasiparticles, and saves the data structure
%   containing the fit, plots the data and the fit, and saves the data 
%   structure containing the fit.

if ~exist('data', 'var')
    data = loadMeasurementData;
end
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable = selectDepDataVars(data, true);
    if isempty(data_variable)
        return
    end
    data_variable = data_variable{1};
end

if exist('T1_expt', 'var')
    T1_set = true;
    T1 = T1_expt;
else
    T1_set = false;
    T1 = [];
end

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);

if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

% Plot 1D data.
if length(dep_rels) == 1
    indep_name = dep_rels{1};
    indep_vals = data.(indep_name);

    f = QPDecayFit(indep_vals(:), dep_vals(:), T1_set, T1);
    
    xunits = getUnits(data, indep_name);
    yunits = getUnits(data, data_variable);
    
    ci = confint(f);
    ae = max([abs(ci(1, 1) - f.a), abs(ci(2, 1) - f.a)]);
    astr = ['a = ', num2str(f.a, 4), ' ± ', num2str(ae, 3), yunits];
    be = max([abs(1 / ci(1, 2) - 1 / f.b), abs(1 / ci(2, 2) - 1 / f.b)]);
    invbstr = ['1/b = ', num2str(1 / f.b, 4), ' ± ', num2str(be, 3),...
        xunits];
    ce = max([abs(ci(1, 3) - f.c), abs(ci(2, 3) - f.c)]);
    cstr = ['c = ', num2str(f.c, 4), ' ± ', num2str(ce, 3), yunits];
    de=max([abs(1 / ci(1, 4) - 1 / f.d), abs(1 / ci(2, 4) - 1 / f.d)]);
    invdstr = ['1/d = ', num2str(1 / f.d, 4), ' ± ', num2str(de, 3),...
        xunits];
    nqpe=max([abs(ci(1, 5) - f.n_qp), abs(ci(2, 5) - f.n_qp)]);
    nqpstr = ['n_{qp} = ', num2str(f.n_qp, 4), ' ± ', num2str(nqpe, 3)];
    full_title = {[strrep(filename, '_', '\_'), ext,...
            ' [', data.Timestamp, ']'],...
            [strrep(data_variable, '_', ' '),' = a * exp(n_{qp}(exp(-d * ',strrep(indep_name, '_', ' '),')-1)-b * ',...
            strrep(indep_name, '_', ' '), ') + c'],...
            [astr, '; ', invbstr, '; ' cstr],...
            [invdstr, '; ',nqpstr]};

    name = ['Fitted_', data_variable];
    data.units.(name) = data.units.(data_variable);
    data.rels.(name) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = name;
    data.plotting.(name).plot_title = full_title;
    data.(name) = f(indep_vals);
    
    % Plot an errorbar graph.
    if isfield(data, 'error') && isfield(data.error, data_variable)
        createFigure('right');
        plotErrorbar(indep_vals, dep_vals, data.error.(data_variable))
        hold on
        plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
        hold off
        legend('data', 'fit')
        xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
        ylabel([strrep(data_variable, '_', ' ') yunits], 'FontSize', 14)
        title(full_title, 'FontSize', 10)
        savePlot(fullfile(plts_path, [filename, '_', data_variable,...
            '_expfit_errorbar']));
    end

    % Plot a simple 1D graph.
    createFigure;
    plotSimple(indep_vals, dep_vals, '.')
    hold on
    plot(indep_vals, f(indep_vals), 'r', 'Linewidth', 2)
    hold off
    legend('data', 'fit')
    xlabel([strrep(indep_name, '_', ' '), xunits], 'FontSize', 14)
    ylabel([strrep(data_variable, '_', ' '), yunits], 'FontSize', 14)
    title(full_title, 'FontSize', 10)
    savePlot(fullfile(plts_path, [filename, '_', data_variable,...
        '_expfit_simple']));

elseif length(dep_rels) == 2 % Plot 2D data.
%    if  ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout'))
    if ~isempty(strfind(dep_rels{1}, 'Time')) ||...
            ~isempty(strfind(dep_rels{1}, 'Duration')) ||...
            ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout')) ||...
            ~isempty(strfind(dep_rels{1}, 'QB_Drive_to_RO'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
%    elseif ~isempty(strfind(dep_rels{2}, 'Qubit_Drive_to_Readout'))
    elseif ~isempty(strfind(dep_rels{2}, 'Time')) ||...
            ~isempty(strfind(dep_rels{2}, 'Duration')) ||...
            ~isempty(strfind(dep_rels{2}, 'Qubit_Drive_to_Readout')) ||...
            ~isempty(strfind(dep_rels{2}, 'QB_Drive_to_RO'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
        flip_fit = false;
    else
        error(['The data does not appear to depenend on ',...
            '''Qubit Drive to Readout'''])
    end

    dep_units = data.units.(data_variable);
    indep2_units = data.units.(indep_name2);

    
    if ~T1_set
        time_const = zeros(length(indep_vals1), 3);
        amplitude = zeros(length(indep_vals1), 3);
        offset = zeros(length(indep_vals1), 3);
        qp_time_const = zeros(length(indep_vals1), 3);
        n_qp_avg = zeros(length(indep_vals1), 3);
        fitted = zeros(size(dep_vals));
        for k = 1:length(indep_vals1)
            f = QPDecayFit(indep_vals2, dep_vals(k, :), T1_set, T1);
            ci = confint(f);
            amplitude(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
            time_const(k, :) = [1/f.b, 1 / (ci(1, 2) - f.b),...
                                   1 / (f.b - ci(2, 2))];
            offset(k, :) = [f.c, f.c - ci(2, 3), ci(1, 3) - f.c];
            qp_time_const(k,:) = [1/f.d, 1 / (ci(1, 4) - f.d),...
                                   1 / (f.d - ci(2, 4))];
            n_qp_avg(k, :) = [f.n_qp, f.n_qp - ci(2, 5), ci(1, 5) - f.n_qp];                   
            fitted(k, :) = f(indep_vals2);
        end
    else
        amplitude = zeros(length(indep_vals1), 3);
        offset = zeros(length(indep_vals1), 3);
        qp_time_const = zeros(length(indep_vals1), 3);
        n_qp_avg = zeros(length(indep_vals1), 3);
        fitted = zeros(size(dep_vals));
        for k = 1:length(indep_vals1)
            f = QPDecayFit(indep_vals2, dep_vals(k, :), T1_set, T1);
            ci = confint(f);
            amplitude(k, :) = [f.a, f.a - ci(2, 1), ci(1, 1) - f.a];
            offset(k, :) = [f.c, f.c - ci(2, 2), ci(1, 2) - f.c];
            qp_time_const(k,:) = [1/f.d, 1 / (ci(1, 3) - f.d),...
                                   1 / (f.d - ci(2, 3))];
            n_qp_avg(k, :) = [f.n_qp, f.n_qp - ci(2, 4), ci(1, 4) - f.n_qp];                   
            fitted(k, :) = f(indep_vals2);
        end
    end
    if flip_fit
        fitted = fitted';
    end
    
    name = 'Extracted_Amplitude';
    data.(name) = amplitude(:, 1);
    data.error.(name) = amplitude(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')

    name = 'Extracted_Offset';
    data.(name) = offset(:, 1);
    data.error.(name) = offset(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    if ~T1_set
        name = 'Extracted_Time_Constant';
        data.(name) = time_const(:, 1);
        data.error.(name) = time_const(:, 2:3);
        data.units.(name) = indep2_units;
        data.rels.(name){1} = indep_name1;
        data.dep{length(data.dep) + 1} = name;
        plotDataVar(data, name, 'errorbar')
    end
    
    name = 'Extracted_QP_Time_Constant';
    data.(name) = qp_time_const(:, 1);
    data.error.(name) = qp_time_const(:, 2:3);
    data.units.(name) = indep2_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    name = 'Extracted_Average_QP_Number';
    data.(name) = n_qp_avg(:, 1);
    data.error.(name) = n_qp_avg(:, 2:3);
    data.units.(name) = dep_units;
    data.rels.(name){1} = indep_name1;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'errorbar')
    
    name = ['Fitted_', data_variable];
    data.(name) = fitted;
    data.units.(name) = dep_units;
    data.rels.(name) = dep_rels;
    data.dep{length(data.dep) + 1} = name;
    plotDataVar(data, name, 'pixelated');
    
    saveMeasData(data, [filename, '_', data_variable, '_expfit'])
end
end

function f = QPDecayFit(x, y, t1_bool, t1_val)
    opts = fitoptions('Method', 'NonlinearLeastSquares',...
                      'Robust', 'LAR',...
                      'MaxFunEvals',10000,...
                      'MaxIter',10000,...
                      'DiffMaxChange',1,...
                      'Algorithm','Trust-Region');
    % opts.Display = 'Off';
    opts.TolX = 1e-9;
    opts.TolFun = 1e-9;
    
    if ~t1_bool
        opts.StartPoint = [y(1) - y(end), 1 / (max(x) - min(x) + 9 * eps),...
                y(end), 0.5 / (max(x) - min(x) + 9 * eps), 1];
        opts.Lower = [min(y), 0, min(y), 0, 0];
        opts.Upper = [max(y), 100*max(x), max(y), 100*max(x), 10];
        f = fit(x(:), y(:), 'a * exp(n_qp * (exp(-d * x) - 1) - b*x) + c', opts);
    else
        opts.StartPoint = [y(1) - y(end), y(end), 0.5 / (max(x) - min(x) + 9 * eps), 1];
        opts.Lower = [min(y), min(y), 0, 0];
        func = ['a * exp(n_qp * (exp(-d * x) - 1) -x/',num2str(t1_val),') + c'];
        f = fit(x(:), y(:), func, opts);
    end
end