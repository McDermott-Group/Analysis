function plotSimple(indep, dep)
%plotSimple Simple 1D data plot.
%
%   plotSimple(INDEP, DEP) plots 1D data. INDEP is an indepedent 
%   variable (x coordinate), DEP is a dependent variable (y coordinate).

plot(indep, dep, '.-', 'LineWidth', 1, 'MarkerSize', 15)

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

axis([xmin xmax ymin ymax])
grid on
set(gca, 'FontSize', 14);

end