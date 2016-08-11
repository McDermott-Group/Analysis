function data = extractKappa
%extractKappa Extract cavity kappa using
% two data sets: Switching Probability vs Readout Amplitude and Switching
% Probability vs Displacement to Fast Pulse.


% Select files for quantum efficiency estimation.
[filenames, pathnames, status] = selectMeasurementDataFile(2,...
    {['Select the first data file containing Switching Probability vs ',...
     'Readout Amplitude...'],...
     ['Select the second data file containing Switching Probability vs',...
     ' Displacement to Readout...']});
if ~status
    return
end

% Read probability data file, convert the variable names, and define
% the units.
file1 = fullfile(pathnames{1}, filenames{1});
file2 = fullfile(pathnames{2}, filenames{2});
data1 = processMeasurementData(importMeasurementData(file1));
data2 = processMeasurementData(importMeasurementData(file2));

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathnames{1});

[~, filename1] = fileparts(filenames{1});
[~, filename2] = fileparts(filenames{2});

sw_prob = 'Switching_Probability';
ro_ampl = 'Readout_Amplitude';
disp2fp = 'Displacement_to_Fast_Pulse';

if ~isfield(data1, sw_prob) || ~isfield(data1, ro_ampl) ||...
        length(data1.(sw_prob)) ~= length(data1.(ro_ampl)) ||...
        length(data1.rels.(sw_prob)) ~= 1 ||...
        ~strcmp(data1.rels.(sw_prob){1}, ro_ampl)
    error(['No Switching Probaility vs Readout Amplitude 1D dataset ',...
        'is found.'])
end

if ~isfield(data2, sw_prob) || ~isfield(data2, disp2fp) ||...
        length(data2.(sw_prob)) ~= length(data2.(disp2fp)) ||...
        length(data2.rels.(sw_prob)) ~= 1 ||...
        ~strcmp(data2.rels.(sw_prob){1}, disp2fp)
    error(['No Switching Probaility vs Displacement to Fase Pulse 1D ',...
        'dataset is found.'])
end

Pkappa = data2.(sw_prob);
ampl = nan(size(Pkappa));
for k = 1:length(Pkappa)
    ampl(k) = Pbright(data1.(ro_ampl), data1.(sw_prob), Pkappa(k));
end

data = data2;
n_photon_name = 'Number_of_Photons';
data.(n_photon_name) = ampl.^2;
data.units.(n_photon_name) = 'Arb. Units';
data.rels.(n_photon_name) = data.rels.(sw_prob);
data.dep{length(data.dep) + 1} = n_photon_name;
data.plotting.(n_photon_name).plot_title =...
    {['Calibration: ', strrep(filenames{1}, '_', '\_'),...
    ' [', data1.Timestamp, ']'],...
    ['Data: ', strrep(filenames{2}, '_', '\_'),...
    ' [', data2.Timestamp, ']']};
data.plotting.(n_photon_name).plot_filename =...
        fullfile(plts_path, [filename1, '_',...
        filename2, '_',  n_photon_name]);
    
plotDataVar(data, n_photon_name);
end

function x0 = Pbright(x, y, y0)
x0 = (Pbright_plus(x, y, y0) + Pbright_minus(x, y, y0)) / 2;
end

function x0 = Pbright_plus(x, y, y0)
    x0 = min(x(y >= y0));
    if isempty(x0)
        x0 = 1;
    end
end

function x0 = Pbright_minus(x, y, y0)
    x0 = max(x(y <= y0));
    if isempty(x0)
        x0 = 0;
    end
end