function plotSemilog(indep, dep)
%plotSemilog Semilog 1D data plot.
%
%   plotSemilog(INDEP, DEP) plots 1D data in a graph with the semilog y
%   axis. INDEP is an indepedent  variable (x coordinate), DEP is
%   a dependent variable (y coordinate).

semilogy(indep, dep, '.-', 'LineWidth', 1, 'MarkerSize', 15)

xmin = min(indep);
xmax = max(indep);
if xmax == xmin
    xmax = Inf;
end

xlim([xmin xmax])
grid on
set(gca, 'FontSize', 14);

end