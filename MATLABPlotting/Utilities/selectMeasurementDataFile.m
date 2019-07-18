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

MAX_NUMBER_OF_FILES = 100;
% status is true if the selection was not interrupted and non-zero number
% of files were chosen.
status = true;

if exist('number_of_files', 'var') && exist('window_titles', 'var') &&...
        number_of_files ~= length(window_titles)
    error(['The number of files to be selected and the length of the ',...
        ' cell that contains the selection window titles do not match.']) 
end

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
if strcmp(pathname, '') || ~exist(pathname, 'dir')
    if exist('Z:\mcdermott-group\Data\Matched JPM Photon Counting', 'dir')
        pathname = 'Z:\mcdermott-group\Data\Matched JPM Photon Counting';
    elseif exist('Z:\Data\Matched JPM Photon Counting', 'dir')
        pathname = 'Z:\Data\Matched JPM Photon Counting';
    end
end

if ~exist('number_of_files', 'var')
    number_of_files = MAX_NUMBER_OF_FILES;
end
filenames = cell(number_of_files, 1);
pathnames = cell(number_of_files, 1);

for k = 1:number_of_files
    % Open the user interface to select a file.
    if ~exist('window_titles', 'var')
        window_title = 'Select a data file...';
    else
        window_title = window_titles{k};
    end
    window_filename = fullfile(pathname, ['MeasurementData_',...
        num2str(k, '%03d')]);
    [filename, pathname] = uigetfile({'*.txt;*.mat;*.hdf5', 'Data Files';
              '*.txt', 'Text Files';
              '*.mat', 'MATLAB Files';
              '*.hdf5', 'HDF5 Files';
              '*.*', 'All Files' }, window_title, window_filename, 'MultiSelect', 'on');

    % Check whether the file selection process should be interrupted.
    if isnumeric(filename)
        if number_of_files < MAX_NUMBER_OF_FILES || k == 1
            status = false;
            return
        end
        break
    else 
        if(iscell(filename) && length(filename) > 1)
            for f = 1:length(filename)
                filenames{k-1+f} = filename{f};
                pathnames{k-1+f} = pathname;
            end
            k = k + length(filename) - 1
        else
            filenames{k} = filename;
            pathnames{k} = pathname;
        end
    end

    % Save the last path used.
    fid = fopen(fullfile(tempdir,...
        'plotMeasurementData_last_pathname.txt'), 'w');
    if fid ~= -1
        fprintf(fid, '%s', pathname);
        fclose(fid);
    end
end

filenames = filenames(~cellfun('isempty', filenames));
pathnames = pathnames(~cellfun('isempty', pathnames));

if ~isempty(filenames)
    if (isfinite(number_of_files) &&...
            length(filenames) == number_of_files) ||...
            ~isfinite(number_of_files)
        status = true;
    end
    if length(filenames) == 1
        filenames = filenames{1};
        pathnames = pathnames{1};
    end
end