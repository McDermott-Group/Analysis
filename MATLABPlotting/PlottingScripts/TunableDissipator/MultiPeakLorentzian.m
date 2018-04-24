function [fitresult, gof] = MultiPeakLorentzian(CutNum, n, InitialPt)
%By Naveen
% Fits multipeak lorentzian data
if CutNum == 0 
    data = loadMeasurementData;
    Freq = data.Cavity_Frequency;
    Occ1 = data.Average_Amplitude;
    Occ1 = Occ1.^2;
    [xData, yData] = prepareCurveData(Freq, Occ1);
else
    data = loadMeasurementData;
    Freq = data.Qubit_Frequency;
    Occ1 = data.Occupation(CutNum,:);
    [xData, yData] = prepareCurveData(Freq, Occ1);
end


if nargin < 1
%     InitialPt = [0.001,0.001,0.001,0.001,0.001,6.632,6.629]; CutNum = 1; n = 2;
    InitialPt = [1e-11,1e-11,1e-11,0.001,0.001,6.632,6.629]; CutNum = 1; n = 2;
    fprintf('FirstCut_doubleLorentzianFit_if2Ddata')
elseif nargin < 2
    InitialPt = [0.001,0.001,0.001,0.001,0.001,4.949,4.946]; CutNum = 1; n = 2;
    fprintf('One Arugement supplied_firstCut_doubleLorentzian')
elseif nargin < 3
     for i = 1:n
        fq(i)=4.949-0.003*(i-1);
     end 
     InitialPt = [ones(1,2*n)*0.001 0.001 fq]
end
   



% Set up fittype and options.
str = '';
for i = 1:n
str = strcat(str,'+a',num2str(i),'*g',num2str(i),'/(2*pi*((x-x',num2str(i),')^2+g',num2str(i),'^2/4))');
end
str = strcat(str,'+c');

ft = fittype( str, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = InitialPt;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
set(h,'Color',[0.2,0.5,0.8]);
set(h,'MarkerSize',10)
legend( h, 'Occupation vs. Qubit Frequency', 'Multipeak Lorentzian Fit', 'Location', 'NorthEast' );

set(h, 'LineWidth',2)
set(h,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);



% Label axes
xlabel Qubit\_Frequency
ylabel Occupation
grid on
if (fitresult.a2*fitresult.a1)>0
    nbar = (fitresult.a1/fitresult.a2-1)^-1
else
    nbar = (-fitresult.a1/fitresult.a2-1)^-1    
end
end


