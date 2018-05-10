function StarkCalibrationDataPlotFromXLS(nbar, GammaPhi)
%CREATEFIT(NBAR,GAMMAPHI)
%
%  Data for 'Linear fit' fit:
%      X Input : nbar
%      Y Output: GammaPhi
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

%  By Naveen, March 28, 2018

a = 5; %Row number start
b = 9; %Row number end
n_col = 10; %Column number for nbar
GammaPhi_col = 8; %Column number for GammaPhi
OFF = 1
if exist('nbar', 'var') && exist('GammaPhi', 'var')
    nbar = nbar;
    GammaPhi = GammaPhi;
elseif exist('OFF','var')
    if OFF == 0
        [~,~,dat]=xlsread('Z:\mcdermott-group\users\Naveen\Projects\Dissipator\EXCEL\April 2018\GammaPhiVsNbars_30thApril2018.xlsx', 'Dissipator ON')
    else
        [~,~,dat]=xlsread('Z:\mcdermott-group\users\Naveen\Projects\Dissipator\EXCEL\April 2018\GammaPhiVsNbars_30thApril2018.xlsx', 'Dissipator OFF')   
    end
    nbar1 = dat(a:b,n_col);
    nbar1_err = dat(a:b,n_col+1);
%     nbar1_err = dat(a:b,n_col+1);
    GammaPhi1 = dat(a:b,GammaPhi_col);
    GammaPhi1_err = dat(a:b,GammaPhi_col+1);
    nbar = cell2mat(nbar1);
    GammaPhi = cell2mat(GammaPhi1);
    nbar_err = cell2mat(nbar1_err);
    GammaPhi_err = cell2mat(GammaPhi1_err);
end



%% Fit: "Linear Fit".
[xData, yData] = prepareCurveData( nbar, GammaPhi );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft )

% Plot fit with data.
figure( 'Name', 'Linear Fit' );
fitresult
h = plot( fitresult, xData, yData );
hold on
% err = errorbar( xData, yData, GammaPhi_err, GammaPhi_err,  nbar_err,  nbar_err)


Slope_confint = confint(fitresult, 0.6827)
GammaPhi_confint = confint(fitresult, 0.6827)

std_Slope_confint = (Slope_confint(2,1) - Slope_confint(1,1))/2;
std_GammaPhi_confint = (GammaPhi_confint(2,2) - GammaPhi_confint(1,2))/2;

% Label axes, title and Equation
if fitresult.p2 < 0
    Equation1 = ['$\Gamma_{\phi E} = ',num2str(round(fitresult.p1,3)),'\hspace{0.1cm} \overline{n} - ',num2str(round(-fitresult.p2,3)),'$'];
    Equation2 = ['$\Delta$','Slope  = ','$ \pm',num2str(std_Slope_confint), '$'];   
    Equation3 = ['$\Delta$','Intercept  = ','$ \pm',num2str(std_GammaPhi_confint), '$'];    
    Eqt5 = text('interpreter','latex','string',{Equation1, Equation2, Equation3});
    set(Eqt5,'Position',[0.302,0.175, 0]);
else
    Equation1 = ['$\Gamma_{\phi E} = ',num2str(round(fitresult.p1,3)),'\hspace{0.1cm} \overline{n} + ',num2str(round(fitresult.p2,3)),'$'];
    Equation2 = ['$\Delta$','Slope  = ','$ \pm',num2str(std_Slope_confint), '$'];   
    Equation3 = ['$\Delta$','Intercept  = ','$ \pm',num2str(std_GammaPhi_confint), '$'];    
    Eqt5 = text('interpreter','latex','string',{Equation1, Equation2, Equation3});
%     set(Eqt5,'Position',[0.122,0.225, 0]);
end 



if OFF == 0
title(['$(\Gamma_{\phi E})_{ON} \hspace{0.3cm} vs \hspace{0.3cm} \overline{n}','$'],'interpreter','latex','FontSize',16)
else
title(['$(\Gamma_{\phi E})_{OFF} \hspace{0.3cm} vs \hspace{0.3cm} \overline{n}','$'],'interpreter','latex','FontSize',16)
end

xlabel(['$\overline{n}','$'],'interpreter','latex','FontSize',16)

ylabel(['$\Gamma_{\phi E} \hspace{0.1cm}  (\mu sec^{-1})','$'],'interpreter','latex','FontSize',16)

h(1).LineStyle = 'none'; 
h(1).Marker = 'o'; 
h(1).MarkerSize = 6;
h(1).MarkerFaceColor = [0, 0.447 , 0.741];
h(1).MarkerEdgeColor = 'none';
h(2).LineWidth = 2;
h(2).Color = 'k';
leg = legend( h, 'Data MIT Procedure', 'Linear Fit', 'Location', 'SouthEast' );
fitresult.p2/fitresult.p1
% err.LineStyle = 'none'; 
grid on
axis tight


