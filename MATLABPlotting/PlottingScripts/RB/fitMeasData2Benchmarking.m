function [out]=fitMeasData2Benchmarking(loadedRBData)

% Select a file if loadedRBData doesn't exist
if ~exist('loadedRBData', 'var')
    data = loadMeasurementData;
else
    data = loadedRBData;
end

% Check loaded data for possibility of generating errorbars
if isfield(data, 'Interleaved_Probabilities')
    errors = 1;
else
    errors = 0;
end

%Find P1 for Interleaved and no interleave for each number of Cliffords
p1_array = zeros(length(data.Number_of_Cliffords), 3);%NoInter, Inter, NoC
p1_array(:,1) = 1 - data.No_Interleaved_Probability;
p1_array(:,2) = 1 - data.Interleaved_Probability;
p1_array(:,3) = data.Number_of_Cliffords;

%Plot
%figure;hold on
%plot(p1_array(:,3), p1_array(:,1), '.-b')
%plot(p1_array(:,3), p1_array(:,2), '.-r')
%legend('No Interleave','Interleave')

% Fit the two curves for depolaring parameter

[out.no_int_x, out.no_int_y]=prepareCurveData(p1_array(:,3), p1_array(:,1));
[out.fr_No_Interleave, ~] = fidelityFit(out.no_int_x, out.no_int_y);

[out.int_x, out.int_y]=prepareCurveData(p1_array(:,3), p1_array(:,2));
[out.fr_Interleave, ~] = fidelityFit(out.int_x, out.int_y);

r_c_est = 0.5*(1-out.fr_Interleave.b/out.fr_No_Interleave.b);
out.fidelity_gate = 1 - r_c_est;

fr_i_conf = confint(out.fr_Interleave);
fr_no_conf = confint(out.fr_No_Interleave);
std_i_b = (fr_i_conf(2,2) - fr_i_conf(1,2))/4;
std_no_b = (fr_no_conf(2,2) - fr_no_conf(1,2))/4;

term1 = std_i_b/out.fr_No_Interleave.b;
term2 = out.fr_Interleave.b/out.fr_No_Interleave.b^2 * std_no_b;
out.sigma_r_c_est = 0.5 * sqrt(term1^2 + term2^2);

% Plotting
createFigure; hold on;

%plot data with or without errorbars based on whether or not the datafile
% contains the required information
if errors
    out.no_int_std = std(data.No_Interleaved_Probabilities,1,2);
    out.int_std = std(data.Interleaved_Probabilities,1,2);
    plotErrorbar(out.no_int_x,out.no_int_y,out.no_int_std);
    plotErrorbar(out.int_x,out.int_y,out.int_std);
else
    plotSimple(out.no_int_x,out.no_int_y);
    plotSimple(out.int_x,out.int_y);
end

%plot the fits on top of the data
powerfunc=@(a,b,c,x) a.*b.^x+c;
m = min(data.Number_of_Cliffords):0.01:max(data.Number_of_Cliffords);
no_int_fit_y = powerfunc(out.fr_No_Interleave.a,out.fr_No_Interleave.b,...
                         out.fr_No_Interleave.c,m);
int_fit_y = powerfunc(out.fr_Interleave.a,out.fr_Interleave.b,...
                         out.fr_Interleave.c,m);

ax=gca;
ax.ColorOrderIndex=1;
plotDots(m,no_int_fit_y);
plotDots(m,int_fit_y);

[~,filename,ext] = fileparts(data.Filename);
filename=strrep(filename,'_','\_');
fid_str = [data.Interleaving_Gate,' Gate Fidelity: ',...
           num2str(out.fidelity_gate),'\pm',num2str(out.sigma_r_c_est)];
title({[filename,ext,': ',data.Timestamp],...
       fid_str});
xlabel('Number of Cliffords')
ylabel('Sequence Fidelity')
leg=legend('Interleaved Gate: None',...
      ['Interleaved Gate: ',data.Interleaving_Gate]);
set(leg,'Location','best')
end

function [fitresult, gof] = fidelityFit(xData, yData)

% Set up fittype and options.
ft = fittype( 'a*b^x + c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.4 0.9 0.5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data
%figure(12);hold on
%plot( fitresult, xData, yData );
end