function extractGamma(data_variable)
% EXTRACTGAMMA Extract JPM switching rate from a set of data files of FPA
% vs FPT.


if ~exist('data_variable', 'var')
    error(['No dependent data variable to fit the exponent to is ',...
        'given as an input argument.'])
end

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

% Read the data files
try
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        raw_data{k} = processMeasurementData(importMeasurementData(file));
        fit_data{k} = fitMeasData2Exponent(data_variable, raw_data{k});
    end
catch
    error('The selected files contain unmatched data.')
end

[legend_entries, ~] = selectLegendEntries(raw_data, data_variable);

createFigure; hold on;
colord = get(gca, 'ColorOrder');
[m,~] = size(colord);
for k=1:length(filenames)
    crow = rem(k,m);
    if crow == 0
        crow = m;
    end
    col = colord(crow,:);
    plot(fit_data{k}.Fast_Pulse_Amplitude, ...
        fit_data{k}.Extracted_Decay_Constant*1e3, '.', 'Color', col, ...
        'MarkerSize', 15);
end
set(gca, 'YScale', 'log');
xlabel('Fast Pulse Amplitude', 'FontSize', 15);
ylabel('Switching Rate (MHz)', 'FontSize', 15);
grid on; set(gca, 'FontSize', 14);
hold off;
legend(legend_entries, 'Interpreter', 'None');

