function deleteDirPlots(pathname, plots_dir_name)
%deleteDirPlots   Delete plot folder in the given path.
%
%   deleteDirPlots(PATHNAME, PLOTS_DIR_NAME) deletes folder PLOTS_DIR_NAME in
%   folder PATHNAME if the folder is empty. PLOTS_DIR_NAME is set
%   to Plots by default.


if ~exist('plots_dir_name', 'var')
    plots_dir_name = 'Plots';
end

plts_path = fullfile(pathname, plots_dir_name);
if exist(plts_path, 'dir')
    rmdir(pathname, plots_dir_name)
end

end

