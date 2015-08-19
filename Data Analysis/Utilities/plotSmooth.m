function plotSmooth(indep1, indep2, dep)
%plotSmooth Smooth 2D data plot.
%
%   plotSmooth(INDEP1, INDEP2, DEP) plots 2D data. INDEP1 and
%   INDEP2 are indepnedent 1D coordinates, x and y correspondingly. DEP is 
%   a dependent variable (z coordinate).

[Ind1, Ind2] = ndgrid(indep1, indep2);
hndl = surf(Ind1, Ind2, dep);
set(gca, 'View', [0 90])
set(hndl, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong');
zmin = min(dep(:));
zmax = max(dep(:));
if zmin == zmax
    zmax = Inf;
end
axis([min(indep1), max(indep1), min(indep2), max(indep2), zmin, zmax, caxis])
colormap(jet)
colorbar

end