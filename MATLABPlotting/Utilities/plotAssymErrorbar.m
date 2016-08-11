function plotAssymErrorbar(indep, dep, l_error, u_error)
%plotAssymErrorbar   Errorbar 1D data plot.
%
%   plotAssymErrorbar(INDEP, DEP, L_ERROR, U_ERROR) plots 1D data. INDEP is
%   an indepedent variable (x coordinate), DEP is a dependent variable
%   (y coordinate). L_ERROR is the length of the lower the error bar and
%   and U_ERROR is the length of the upper error bar.


errorbar(indep, dep, l_error, u_error, '.',...
    'LineWidth', 1, 'MarkerSize', 15)

axis tight
grid on
set(gca, 'FontSize', 14)

end