function plotDots(indep, dep)
%plotDots Dotted 1D data plot.
%
%   plotDots(INDEP, DEP) plots 1D data. INDEP is an indepedent 
%   variable (x coordinate), DEP is a dependent variable (y coordinate).

plot(indep, dep, '+', 'MarkerSize', 5)

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