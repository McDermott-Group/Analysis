function savePlot(full_filename)
%savePlot   Save current figure plot in a png file.
%
%   savePlot(FULL_FILENAME) saves the plot in a current figure to a .png
%   file. FULL_FILENAME specifies the name of the png file and its location.

try
    saveas(gca, full_filename, 'png')
catch
    disp(['Plot ''', full_filename, ''' could not be saved. ',...
            'The automatically created file name is probably too long.'])
end

end