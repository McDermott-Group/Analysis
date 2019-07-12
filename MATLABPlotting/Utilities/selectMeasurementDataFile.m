function [filenames, pathnames, status] =...
        selectMeasurementDataFile(number_of_files, window_titles)
%SelectMeasurementDataFile  Open the file selection dialog to select a set
%of files. 
%
%   FILENAMES, PATHNAMES, STATUS = selectMeasurementDataFile(NUMBER_OF_FILES,
%   WINDOW_TITLES) opens the file selection dialog NUMBER_OF_FILES times. 
%   WINDOW_TITLES should be a cell of strings. Each string specifies
%   the window title of the n-th file selection dialog window.
%   The length of WINDOW_TITLES should be equal to NUMBER_OF_FILES.
%   The function returns cell FILENAMES contaning the chosen filenames, 
%   cell PATHNAMES containg corresponding pathes and STATUS flag that is
%   true if the file selection process is finised as expected.
%
%   If NUMBER_OF_FILES is not specified an abitrary number of files could
%   be selected. The selection process will stop upon pressing [Cancel]
%   button. If FILENAMES and PATHNAMES are cells containg only one string
%   each they will be converted to a simple string.


% Retrive the path that was used last time.
pathname = '';
fid = fopen(fullfile(tempdir, 'plotMeasurementData_last_pathname.txt'),...
    'r');
if fid ~= -1
    pathname = fgetl(fid);
    fclose(fid);
end

% If the path hasn't be retrieved or does not exist, predefine it with
% one of the potentially existing path values.
% if strcmp(pathname, '') || ~exist(pathname, 'dir')
%     if exist('Z:\mcdermott-group\Data\Matched JPM Photon Counting', 'dir')
%         pathname = 'Z:\mcdermott-group\Data\Matched JPM Photon Counting';
%     elseif exist('Z:\Data\Matched JPM Photon Counting', 'dir')
%         pathname = 'Z:\Data\Matched JPM Photon Counting';
%     end
% end

% Open the user interface to select a file.
window_title = 'Select a data file...';
window_filename = fullfile(pathname, ['MeasurementData_',...
    num2str(1, '%03d')]);
[filenames, pathnames] = uigetfile({'*.txt;*.mat;*.hdf5', 'Data Files';
          '*.txt', 'Text Files';
          '*.mat', 'MATLAB Files';
          '*.hdf5', 'HDF5 Files';
          '*.*', 'All Files' }, window_title, window_filename, 'MultiSelect', 'on');

% Save the last path used.
fid = fopen(fullfile(tempdir,...
    'plotMeasurementData_last_pathname.txt'), 'w');
if fid ~= -1
    fprintf(fid, '%s', pathname);
    fclose(fid);
end