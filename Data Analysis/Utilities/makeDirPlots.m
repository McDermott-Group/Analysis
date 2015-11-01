function plts_path = makeDirPlots(pathname, plts_dir_name)
%makeDirPlots Create folder Plots in the given path.
%
%   makeDirPlots(PATHNAME, PLTS_DIR_NAME) creates folder PLTS_DIR_NAME in
%   folder PATHNAME if the former does not exist.

if ~exist('plots_dir_name', 'var')
    plts_dir_name = 'Plots';
end

plts_path = fullfile(pathname, plts_dir_name);
if ~exist(plts_path, 'dir')
    mkdir(pathname, plts_dir_name)
end

end

