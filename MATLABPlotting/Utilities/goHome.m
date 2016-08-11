function goHome
%goHome Jump into the .\Git Repositories\MATLAB\Data Analysis folder.

pkg_path = mfilename('fullpath');
seps = strfind(pkg_path, filesep);
pkg_path = pkg_path(1:seps(end-1) - 1);

cd(pkg_path)