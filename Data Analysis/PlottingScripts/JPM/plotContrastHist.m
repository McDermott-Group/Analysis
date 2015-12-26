function plotContrastHist
%plotContrastHist Plot a contrast histogram.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

[pathname, filename, ext] = fileparts(data.Filename);
plts_path = makeDirPlots(pathname);
plot_title = [filename, ext, ' [', data.Timestamp, ']'];

if ~isfield(data, 'Contrast')
    error('''Contrast'' data not found in the selected file.')
end

max_contrast = ceil(100 * max(data.Contrast));

bins = linspace(.5, max_contrast - .5, max_contrast);

createFigure;
hist(100 * data.Contrast, bins)

xlabel('Contrast (%)', 'FontSize', 14);
ylabel('Statistical Frequency', 'FontSize', 14);
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)

savePlot(fullfile(plts_path, [filename, '_hist']));

[~, idx] = max(data.Contrast);
entries = fieldnames(data);
len = length(data.Contrast);
for k = 1:length(entries)
    name = entries{k};
    vals = data.(name);
    if isnumeric(vals) && length(vals) == len
        units = getUnits(data, name);
        if ~isempty(units)
            units = [' ', units(3:end-1)];
        end
        disp([strrep(name, '_', ' '), ' = ', num2str(vals(idx), 6),...
            units, ' (at maximum), ', num2str(vals(end), 6), units,...
            ' (last)'])
    end
end