function [Rn, Ic] = RdiodeToRnormal(Rmeasured)
% [Rn, Ic] = RdiodeToRnormal(Rmeasured)
% Input is the resistance(in kOhms) with diode mode ON on probe station
% This function maps from those resistances to the actual resistances
% measured. And, it outputs actual resistances and corresponding critical
% currents for Aluminium junctions. 

rDiode1 = [1,4.975,9.949,14.770,32.574,19.320,26.064,31.000,33.59,2.314,7.685,13.686,16.069,22.423,23.304,24.64,30.037,32.668];
rNormal1 = [1,4.984,10.003,68.397,14.984,20.03,30.237,50.12,99.418,2.316,7.71,13.843,16.367,23.975,25.294,27.439,43.26,69.622];
%Discarded data points (32.81,74.851)  
rDiode = sort(rDiode1);
rNormal = sort(rNormal1);

Rn = spline(rDiode,rNormal,Rmeasured);
Ic = pi*3.4e-4./(4*Rn*1000)*1e9;
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
ft1 = 'splineinterp';

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
