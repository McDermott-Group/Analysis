function computeXCovariance
%computeXCovariance computes cross-covariance for two JPM measurements.

% Retrive the path that was used last time.
pathname = '';
fid = fopen(fullfile(tempdir, 'plotMeasurementData_last_pathname.txt'), 'r');
if fid ~= -1
    pathname = fgetl(fid);
    fclose(fid);
end

% If the path hasn't be retrieved or does not exist, predefine it with
% one of the potentially existing values.
if strcmp(pathname, '') || ~exist(pathname, 'dir')
    if exist('Z:\mcdermott-group\Data\Matched JPM Photon Counting', 'dir')
        pathname = 'Z:\mcdermott-group\Data\Matched JPM Photon Counting';
    elseif exist('Z:\Data\Matched JPM Photon Counting', 'dir')
        pathname = 'Z:\Data\Matched JPM Photon Counting';
    end
end

% Open the user interface to select a file.
[filename, pathname] = uigetfile({'*.txt', 'Text Files';...
          '*.*', 'All Files' },...
          'Select a data file...', fullfile(pathname, 'MeasurementData_000'));

if isnumeric(filename)
    return
end

% Save the current path.
fid = fopen(fullfile(tempdir, 'plotMeasurementData_last_pathname.txt'), 'w');
if fid ~= -1
    fprintf(fid, '%s', pathname);
    fclose(fid);
end

% Read probability data file, convert the variable names, and define
% the units.
data = processMeasurementData(importMeasurementData(fullfile(pathname, filename)));

% Sanity checks.
if isempty(data)
    error('The selected file does not contain any data.');
elseif ~isfield(data, 'indep') || isempty(data.indep)
    error('The independent (sweep) variables are not specified.');
elseif ~isfield(data, 'dep') || isempty(data.dep)
    error('The dependent (data) variables are not specified.');
elseif ~isfield(data, 'rels') || isempty(data.rels)
    error('The relationships between the dependent (data) and independent (sweep) variables are not specified.');  
else
    for k = 1:length(data.indep)
        if isempty(data.(data.indep{k}))
            error('The selected file does not specify the independent (sweep) variables.');
        end
    end
    for k = 1:length(data.dep)
        if isempty(data.(data.dep{k}))
            error('The selected file does not contain any actual data.');
        end
    end
    RelationshipsFlag = false; 
    for k = 1:length(data.dep)
        if ~isempty(data.(data.dep{k}))
            RelationshipsFlag = true;
        end
    end
    if ~RelationshipsFlag
        error('The relationships between the dependent (data) and independent (sweep) variables are not specified.'); 
    end
end

for data_index = 1:length(data.dep)
    dep_name = data.dep{data_index};
    dep_vals = data.(dep_name);
    dep_distr = data.distr.(dep_name);
    
    if strcmp(dep_distr, 'binomial')
        data.error.(dep_name) = sqrt(dep_vals .* (1 - dep_vals) / data.Number_of_Measurements);
    elseif strcmp(dep_distr, 'normal')
        if isfield(data, [dep_name, '_Std_Dev'])
            data.error.(dep_name) = data.([dep_name, '_Std_Dev']) / sqrt(data.Number_of_Measurements - 1);
        end
    elseif strcmp(dep_distr, 'std')
        continue
    end
end

if ~isfield(data, 'JPM_A_Switching_Probability') || length(data.rels.JPM_A_Switching_Probability) ~= 1 ||...
   ~isfield(data, 'JPM_B_Switching_Probability') || length(data.rels.JPM_B_Switching_Probability) ~= 1 ||...
   ~isfield(data, 'P00') || length(data.rels.P00) ~= 1 ||...
   ~isfield(data, 'P01') || length(data.rels.P01) ~= 1 ||...
   ~isfield(data, 'P10') || length(data.rels.P10) ~= 1 ||...
   ~isfield(data, 'P11') || length(data.rels.P11) ~= 1
    error(['The data should have the following 1D data variables ''JPM_A_Switching_Probability'', ',...
           '''JPM_B_Switching_Probability'', ', '''P00'', ', '''P01'', ', '''P10'', and ', '''P11''',...
           ' and these data variables should depend on an single independent (sweep) variable.'])
end

if ~strcmp(data.rels.P00{1}, data.rels.JPM_A_Switching_Probability{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.JPM_B_Switching_Probability{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P01{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P10{1}) ||...
   ~strcmp(data.rels.P00{1}, data.rels.P11{1})
    error('The independent (sweep) variable should be the same for all probability data variables.')
end

if length(data.rels.P00) > 1
    error('The can only handle 1D data sets.')
end

indep_vals = data.(data.rels.P00{1});
indep_name = data.rels.P00{1};
n = length(indep_vals);
half_n = floor(n/2);

disp('Auto:')
disp(['Cross-Cov(P00, P00)[tB - tA < 0] = ', num2str(xcov(data.P00(1:half_n), data.P00(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P00, P00)[tB - tA > 0] = ', num2str(xcov(data.P00(half_n+2:end), data.P00(half_n+2:end), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P11, P11)[tB - tA < 0] = ', num2str(xcov(data.P11(1:half_n), data.P11(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P11, P11)[tB - tA > 0] = ', num2str(xcov(data.P11(half_n+2:end), data.P11(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('P00 vs P11:')
disp(['Cross-Cov(P00, P11)[tB - tA < 0] = ', num2str(xcov(data.P00(1:half_n), data.P11(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P00, P11)[tB - tA > 0] = ', num2str(xcov(data.P00(half_n+2:end), data.P11(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('P01 vs P00:')
disp(['Cross-Cov(P01, P00)[tB - tA < 0] = ', num2str(xcov(data.P01(1:half_n), data.P00(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P01, P00)[tB - tA > 0] = ', num2str(xcov(data.P01(half_n+2:end), data.P00(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('P10 vs P00:')
disp(['Cross-Cov(P10, P00)[tB - tA < 0] = ', num2str(xcov(data.P10(1:half_n), data.P00(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P10, P00)[tB - tA > 0] = ', num2str(xcov(data.P10(half_n+2:end), data.P00(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('P01 vs P11:')
disp(['Cross-Cov(P01, P11)[tB - tA < 0] = ', num2str(xcov(data.P01(1:half_n), data.P11(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P01, P11)[tB - tA > 0] = ', num2str(xcov(data.P01(half_n+2:end), data.P11(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('P10 vs P11:')
disp(['Cross-Cov(P10, P11)[tB - tA < 0] = ', num2str(xcov(data.P10(1:half_n), data.P11(1:half_n), 0) / (half_n - 1)), '.']);
disp(['Cross-Cov(P10, P11)[tB - tA > 0] = ', num2str(xcov(data.P10(half_n+2:end), data.P11(half_n+2:end), 0) / (half_n - 1)), '.']);

disp('Auto:')
disp('Corr-Coeff(P00, P00)[tB - tA < 0] = ')
disp(corrcoef(data.P00(1:half_n), data.P00(1:half_n)));
disp('Corr-Coeff(P00, P00)[tB - tA > 0] = ')
disp(corrcoef(data.P00(half_n+2:end), data.P00(half_n+2:end)));
disp('Corr-Coeff(P11, P11)[tB - tA < 0] = ')
disp(corrcoef(data.P11(1:half_n), data.P11(1:half_n)));
disp('Corr-Coeff(P11, P11)[tB - tA > 0] = ')
disp(corrcoef(data.P11(half_n+2:end), data.P11(half_n+2:end)));

disp('P00 vs P11:')
disp('Corr-Coeff(P00, P11)[tB - tA < 0] = ')
disp(corrcoef(data.P00(1:half_n), data.P11(1:half_n)));
disp('Corr-Coeff(P00, P11)[tB - tA > 0] = ')
disp(corrcoef(data.P00(half_n+2:end), data.P11(half_n+2:end)));

disp('P01 vs P00:')
disp('Corr-Coeff(P01, P00)[tB - tA < 0] = ')
disp(corrcoef(data.P01(1:half_n), data.P00(1:half_n)));
disp('Corr-Coeff(P01, P00)[tB - tA > 0] = ')
disp(corrcoef(data.P01(half_n+2:end), data.P00(half_n+2:end)));

disp('P10 vs P00:')
disp('Corr-Coeff(P10, P00)[tB - tA < 0] = ')
disp(corrcoef(data.P10(1:half_n), data.P00(1:half_n)));
disp('Corr-Coeff(P10, P00)[tB - tA > 0] = ')
disp(corrcoef(data.P10(half_n+2:end), data.P00(half_n+2:end)));

disp('P01 vs P11:')
disp('Corr-Coeff(P01, P11)[tB - tA < 0] = ')
disp(corrcoef(data.P01(1:half_n), data.P11(1:half_n)));
disp('Corr-Coeff(P01, P11)[tB - tA > 0] = ')
disp(corrcoef(data.P01(half_n+2:end), data.P11(half_n+2:end)));

disp('P10 vs P11:')
disp('Corr-Coeff(P10, P11)[tB - tA < 0] = ')
disp(corrcoef(data.P10(1:half_n), data.P11(1:half_n)));
disp('Corr-Coeff(P10, P11)[tB - tA > 0] = ')
disp(corrcoef(data.P10(half_n+2:end), data.P11(half_n+2:end)));

[P10vsP11, ~] = xcorr(data.P10 - mean(data.P10), data.P11 - mean(data.P11), half_n, 'coeff');
[P01vsP11, lags] = xcorr(data.P01 - mean(data.P01), data.P11 - mean(data.P11), half_n, 'coeff');

figure
plot(lags, P10vsP11, 'r.', lags, P01vsP11, 'b.')
grid on
legend('Corr-Coeff(P10, P11)', 'Corr-Coeff(P01, P11)')

if isfield(data.units, data.rels.P00{1}) && ~isempty(data.units.(data.rels.P00{1}))
    xunits = [' (', data.units.(data.rels.P00{1}), ')'];
else
    xunits = '';
end

xlabel(['Time Lag', xunits], 'FontSize', 14);
ylabel('Cross-Correlation Coefficient', 'FontSize', 14);
title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)

figure
plot(lags, P10vsP11 - P01vsP11, 'b.')
grid on
legend('Corr-Coeff(P10, P11) - Corr-Coeff(P01, P11)')

xlabel(['Time Lag', xunits], 'FontSize', 14);
ylabel('\Delta Cross-Correlation Coefficient', 'FontSize', 14);
title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)

figure
[PAvsPB, lags] = xcorr(data.JPM_A_Switching_Probability - mean(data.JPM_A_Switching_Probability),...
                       data.JPM_B_Switching_Probability - mean(data.JPM_B_Switching_Probability),...
                       half_n, 'coeff');
plot(lags, PAvsPB, 'b.')
grid on
legend('Corr-Coeff(PA, PB)')

xlabel(['Time Lag', xunits], 'FontSize', 14);
ylabel('Cross-Correlation Coefficient', 'FontSize', 14);
title({[filename, ' [', data.Timestamp, ']']}, 'Interpreter', 'none', 'FontSize', 10)

% Show a message box with the experiment parameters.
fields = fieldnames(data);
params = cell(length(fields), 1);
for k = 1:length(fields)
    if isnumeric(data.(fields{k}))
        if length(data.(fields{k})) == 1
            if isfield(data.units, fields{k}) && ~isempty(data.units.(fields{k}))
                units = [' (', data.units.(fields{k}), ')'];
            else
                units = '';
            end
            params{k} = [strrep(fields{k},  '_', ' '), ' = ',...
                        num2str(data.(fields{k})), ' ', units];
        end
    end
    if strcmp(fields{k}, 'Timestamp')
        params{k} = data.(fields{k});
    end
    if strcmp(fields{k}, 'Comments')
        params{k} = ['Comments: ', [data.(fields{k}){:}]];
    end
end
temp_struct.Interpreter = 'tex';
temp_struct.WindowStyle = 'non-modal';
msgbox(params(~cellfun('isempty', params)), 'Parameters', temp_struct)