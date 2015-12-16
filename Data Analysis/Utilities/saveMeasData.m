function saveMeasData(data, filename)
%saveMeasData  Save data in .mat file.
%
%   saveMeasData(DATA, FILENAME) save the DATA structure in the
%   FILENAME file.

% Retrive the path that was used last time.
pathname = '';
fid = fopen(fullfile(tempdir, 'plotMeasurementData_last_pathname.txt'), 'r');
if fid ~= -1
    pathname = fgetl(fid);
    fclose(fid);
end

% If the path hasn't be retrieved or does not exist, predefine it with
% one of the potentially existing path values.
if strcmp(pathname, '') || ~exist(pathname, 'dir')
    if exist('Z:\mcdermott-group\Data\Data Analysis', 'dir')
        pathname = 'Z:\mcdermott-group\Data\Data Analysis';
    elseif exist('Z:\Data\Data Analysis', 'dir')
        pathname = 'Z:\Data\Data Analysis';
    elseif exist('\\afs\physics.wisc.edu\mcdermott-group\Data\Data Analysis', 'dir')
        pathname = '\\afs\physics.wisc.edu\mcdermott-group\Data\Data Analysis';
    end
end

file = fullfile(pathname, filename);

save(file, 'data');