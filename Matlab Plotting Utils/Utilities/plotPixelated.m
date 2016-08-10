function plotPixelated(indep1, indep2, dep)
%plotPixelated  Pixelated 2D data plot.
%
%   plotPixelated(INDEP1, INDEP2, DEP) plots 2D data as a pixelated image.
%   INDEP1 and INDEP2 are independent 1D coordinates, x and y correspondigly.
%   DEP is a dependent variable (z coordinate).

imagesc(indep1, indep2, dep);

set(gca, 'YDir', 'normal')
colormap(jet)
colorbar

end