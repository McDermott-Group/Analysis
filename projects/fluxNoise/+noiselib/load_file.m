function [data] = load_file(pathname, filename)
% Read the data file, convert the variable names, and specify the units.

try
    data = cell(0);
    file = fullfile(pathname, filename);
    data = loadMeasurementData(file);
catch
    error(['The file ' filename ' contains unmatched data.'])
end

end