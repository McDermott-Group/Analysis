function plotPolar(radius, phase, dep)
%plotPolar  Polar 2D data plot.
%
%   plotPolar(RADIUS, PHASE, DEP) plots 2D data in polar coordinates. 
%   RADIUS and PHASE are polar and phase coordinates, correspondingly. 
%   DEP is a dependent variable (z coordinate).

[Radius, Phase] = ndgrid(radius, phase);
[X, Y] = pol2cart(Phase, Radius);

hndl = surf(X, Y, dep);
set(gca, 'View', [0 90])
set(hndl, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong');
zmin = min(dep(:));
zmax = max(dep(:));
if zmin == zmax
    zmax = Inf;
end
axis([-max(radius), max(radius), -max(radius), max(radius), zmin, zmax, caxis])
axis equal
axis off
colormap(jet)
colorbar

end