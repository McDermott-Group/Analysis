function plotPowerSensitivity(data_variable, data)
%plotPowerSensitivity(DATA_VARIABLE, DATA) Plot power sensitivity defined
%as switching probability derivative over number of photons.
%   plotPowerSensitivity(DATA_VARIABLE, DATA) Plot power sensitivity 
%   for DATA_VARIABLE that is normally expected to be switching
%   probability. DATA is optional data structure.

if ~exist('data', 'var')
    % Select a file.
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

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
makeDirPlots(pathname);

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
    
    if isempty(strfind(indep_name, 'Attenuation')) &&...
            isempty(strfind(indep_name, 'Amplitude'))
        error(['The data does not appear to depenend on any ',...
            '''Attenuaation'' or ''Amplitude'' variable.'])
    end
    
    dep_vals = diff(dep_vals) / median(diff(indep_vals));
    dep_vals(end+1) = dep_vals;
        
    pos = strfind(indep_name, 'Attenuation');
    if ~isempty(pos)
        dep_vals = -dep_vals ./ (10.^(-indep_vals / 10) * log(10) / 10);
    else
        pos = strfind(indep_name, 'Attenuation');
        dep_vals = dep_vals ./ indep_vals;
    end
    
    prefix = indep_name(1:pos(1)-2);
    N = photonFluence(data.([prefix, '_Power']), 60,...
        data.([prefix, '_Frequency']), data.Fast_Pulse_Time);

elseif length(dep_rels) == 2 % Plot 2D data.
    mode = true;
    if ~isempty(strfind(dep_rels{1}, 'Attenuation'))
        pos = strfind(dep_rels{1}, 'Attenuation');
        prefix = dep_rels{1}(1:pos(1)-2);
        indep_vals2 = data.(dep_rels{1});
        dep_rels = {dep_rels{2}, dep_rels{1}};
        dep_vals = dep_vals';
    elseif ~isempty(strfind(dep_rels{2}, 'Attenuation'))
        pos = strfind(dep_rels{2}, 'Attenuation');
        prefix = dep_rels{2}(1:pos(1)-2);
        indep_vals2 = data.(dep_rels{2});
    elseif ~isempty(strfind(dep_rels{1}, 'Amplitude')) &&...
            isempty(strfind(dep_rels{1}, 'Fast'))
        prefix = 'Readout';
        mode = false;
        indep_vals2 = data.(dep_rels{1});
        dep_rels = {dep_rels{2}, dep_rels{1}};
        dep_vals = dep_vals';
    elseif ~isempty(strfind(dep_rels{2}, 'Amplitude')) &&...
            isempty(strfind(dep_rels{2}, 'Fast'))
        prefix = 'Readout';
        indep_vals2 = data.(dep_rels{2});
        mode = false;
    else
        error(['The data does not appear to depenend on any ',...
            '''Attenuaation'' or ''Amplitude'' variable.'])
    end

    dep_vals = conv2(1, ones(1, 2) / 2, dep_vals, 'same');
    dep_vals = diff(dep_vals, 1, 2) / median(diff(indep_vals2));
    dep_vals(:, end+1) = dep_vals(:, end);
    if mode
        attn = ones(size(dep_vals, 1), 1) * indep_vals2(:)';
        dep_vals = -dep_vals ./ (10.^(-attn / 10) * log(10) / 10);
    else
        ampl = ones(size(dep_vals, 1), 1) * indep_vals2(:)';
        dep_vals = dep_vals ./ ampl;
    end
    
    N = photonFluence(data.([prefix, '_Power']), 60,...
        data.([prefix, '_Frequency']), data.Fast_Pulse_Time);

    dep_vals(dep_vals < 0) = 0;
    threshold = 20 * median(dep_vals(:));
    dep_vals(dep_vals > threshold) = threshold;
end

full_title = {'Power Sensitivity',...
    [strrep(filename, '_', '\_'), ext, ' [', data.Timestamp, ']']};

name = 'Power_Sensitivity';
data.(name) = dep_vals / N;
data.units.(name) = '';
data.rels.(name) = dep_rels;
data.dep{length(data.dep) + 1} = name;
data.plotting.(name).plot_title = full_title;
plotDataVar(data, name)
end