function [Rn, Ic] = RdiodeToRnormal(Rmeasured)
% [Rn, Ic] = RdiodeToRnormal(Rmeasured)
% Input is the resistance(in kOhms) with diode mode ON on probe station
% This function maps from those resistances to the actual resistances
% measured. And, it outputs actual resistances and corresponding critical
% currents for Aluminium junctions.

rDiode = [1,4.975,9.949,14.770,19.320,26.064,31.000,32.81,33.59];
rNormal = [1,4.984,10.003,14.984,20.03,30.237,50.12,74.851,99.418];
Rn = spline(rDiode,rNormal,Rmeasured);
Ic = pi*3.4e-4./(2*Rn*1000)*1e9;
Rn = {Rn, 'kOhm'};
Ic = {Ic, 'nA'};
valsRn = [''];
valsIc = [''];

for i = 1:length(Rmeasured)
    valsRn = strcat([valsRn,'  ',num2str(Rn{1}(i)),' ',Rn{2},',']);
    valsIc = strcat([valsIc,'  ',num2str(Ic{1}(i)),' ',Ic{2},',']);
end
valsRn
valsIc


%% Fit resistance data
[xData, yData] = prepareCurveData( rDiode, rNormal );

% Set up fittype and options.
ft = 'splineinterp';

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, 'Normalize', 'on' );

% Plot fit with data.
figure( 'Name', 'Fit' );
h = plot( fitresult, xData, yData );
legend( h, 'R_{normal} vs. R_{diode}', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel('R_{Diode} (k\Omega)')
ylabel('R_{Normal} (k\Omega)')
grid on
