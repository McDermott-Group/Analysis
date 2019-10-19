function [fitresult, gof] = fit_lorenzian(x, y)
% Set up fittype and options.
ft = fittype( 'a/(4*pi^2*x.^2+gamma^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [380 1056];

% Fit model to data.
[fitresult, gof] = fit( x, y, ft, opts );
end