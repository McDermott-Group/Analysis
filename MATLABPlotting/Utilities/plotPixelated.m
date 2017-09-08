function plotPixelated(indep1, indep2, dep)
%plotPixelated  Pixelated 2D data plot.
%
%   plotPixelated(INDEP1, INDEP2, DEP) plots 2D data as a pixelated image.
%   INDEP1 and INDEP2 are independent 1D coordinates,
%   x and y correspondigly. DEP is a dependent variable (z coordinate).

diff1 = abs(diff(indep1));
diff2 = abs(diff(indep2));

if min(diff1) < .95 * max(diff1) || min(diff2) < .95 * max(diff2)
    disp(['The pixelated plot could be incorrect due to the ',...
          'ununiform data point spacing.'])
%     close(gcf);
%     return
end

imagesc(indep1, indep2, dep);

set(gca, 'YDir', 'normal')
colormap(jet)
colorbar

end