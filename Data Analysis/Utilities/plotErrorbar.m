function plotErrorbar(indep, dep, dep_error)
%plotErrorbar Errorbar 1D data plot.
%
%   plotErrorbar(INDEP, DEP, DEP_ERROR) plots 1D data. INDEP is an indepedent 
%   variable (x coordinate), DEP is a dependent variable (y coordinate).
%   DEP_ERROR is the error on DEP variable.


errorbar(indep, dep, 1.96 * dep_error, '.', 'LineWidth', 1, 'MarkerSize', 15)

xmin = min(indep);
xmax = max(indep);
if xmax == xmin
    xmax = Inf;
end

ymin = min(dep - 1.96 * dep_error);
ymax = max(dep + 1.96 * dep_error);
if ymin == ymax
    ymax = Inf;
end

if isfinite(ymin) && isfinite(ymax)
    axis([xmin xmax ymin ymax])
else
    xlim([xmin xmax])
end
grid on
set(gca, 'FontSize', 14)

end