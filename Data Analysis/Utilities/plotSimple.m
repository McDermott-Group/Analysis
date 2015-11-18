function plotSimple(indep, dep, style)
%plotSimple Simple 1D data plot.
%
%   plotSimple(INDEP, DEP, STYLE) plots 1D data. INDEP is an indepedent 
%   variable (x coordinate), DEP is a dependent variable (y coordinate),
%   STYLE is the desired line style.

if ~exist('style', 'var')
    style = '.-';
end

plot(indep, dep, style, 'LineWidth', 1, 'MarkerSize', 15)

xmin = min(indep);
xmax = max(indep);
if xmax == xmin
    xmax = Inf;
end

ymin = min(dep);
ymax = max(dep);
if ymin == ymax
    ymax = Inf;
end

if isfinite(ymin) && isfinite(ymax)
    axis([xmin xmax ymin ymax])
else
    xlim([xmin xmax])
end
grid on
set(gca, 'FontSize', 14);

end