Plot_QP_TunnelingGap('Nb_GND_Dev01\', 'Leiden_2020Feb\LIU\',  'Q1\', '03-17-20\', 0:4);
%Plot_QP_TunnelingGap('Nb_GND_Dev01\', 'Leiden_2020Feb\LIU\',  'Q1\', '03-16-20\', 0);

% plotFit(1/0.28e-3, 2e-5);
% plotFit(1/13e-3, 9e-4);
% plotFit(1/1e-3, 5e-5);
%plotFit(1/1.5e-3, 4e-4);
%plotFit(1/1.6e-2, 4e-4);
%plotFit(1/1.8e-2, 0.8);

%function [] = plotFit(gamma, A0)
%    f = 10.^[-1:0.01:3.5]; 
%    A = gamma/4*A0;
%    theoryplot = gamma*A./(4*pi^2*f.^2+gamma^2/4)+();
%    figure(111); hold on; 
%    plot(f,theoryplot, 'DisplayName', ['fit ',num2str(1e3/gamma),'ms ',num2str(A0),'Fidelity'])
%end

function [] = plotFit(gamma, F)
    f = 10.^[-1:0.01:3.5]; 
    t = 50e-6
    theoryplot = F*F*gamma./(gamma*gamma+pi^2*f.^2)+(1-F*F)*t;
    figure(111); hold on; 
    plot(f,theoryplot, 'DisplayName', ['fit ',num2str(1e3/gamma),'ms ',num2str(F),'Fidelity'])
end