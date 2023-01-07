function [fitresult, gof] = fit_psd(f,psd)

[xData, yData] = prepareCurveData(f, psd');

ft = fittype( 'a + b*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-7 -1.5];

[fitresult, gof] = fit( xData, yData, ft, opts );

% figure;plot( fitresult, xData, yData );

end