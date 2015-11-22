function plotDots(indep, dep)
%plotDots   Dotted 1D data plot.
%
%   plotDots(INDEP, DEP) plots 1D data. INDEP is an indepedent 
%   variable (x coordinate), DEP is a dependent variable (y coordinate).

plot(indep, dep, '+', 'MarkerSize', 5)

axis tight
grid on
set(gca, 'FontSize', 14);

end