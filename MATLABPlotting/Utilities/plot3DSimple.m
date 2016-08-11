function plot3DSimple(indep1, indep2, dep, style)
%plot3DSimple Simple 3D scatter plot.
%
%   plotSimple(INDEP1, INDEP2, DEP, STYLE) scatter plot of 3D data. INDEP1 
%   and INDEP2 are the is an indepedent variables (x and y coordinates), 
%   DEP is a dependent variable (y coordinate), STYLE is the desired 
%   line style.

if ~exist('style', 'var')
    style = '-.';
end

[X,Y] = meshgrid(indep1, indep2)

plot3(X, Y, dep, style, 'MarkerSize', 15)
set(gca, 'View', [90, 0])
axis tight
grid on

end