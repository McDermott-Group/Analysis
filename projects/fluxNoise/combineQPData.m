% Analyze_QP_Tunneling('Circ1\', 'Q2\', '07-21-19\', 21, 51);
% Analyze_QP_Tunneling('Circ1\', 'Q3\', '07-26-19\', 0, 210);
% Analyze_QP_Tunneling('Circ1\', 'Q4\', '07-21-19\', 2, 38);
% Analyze_QP_Tunneling('Circ2\', 'Q1\', '07-26-19\', 0, 70);
% Analyze_QP_Tunneling('Circ2\', 'Q3\', '07-25-19\', 0, 691);
% Analyze_QP_Tunneling('Circ2\', 'Q4\', '07-19-19\', 0, 30);

% Analyze_QP_Tunneling('Circ1\', 'Q1\', '08-22-19\', 0, 218);
% Analyze_QP_Tunneling('Circ1\', 'Q1\', '08-29-19\', 0, 99);
% Analyze_QP_Tunneling('Circ1\', 'Q2\', '08-22-19\', 0, 118);
% Analyze_QP_Tunneling('Circ1\', 'Q3\', '08-21-19\', 0, 88);
% Analyze_QP_Tunneling('Circ1\', 'Q4\', '08-22-19\', 0, 204);
% Analyze_QP_Tunneling('Circ1\', 'Q4\', '08-29-19\', 0, 99);
% Analyze_QP_Tunneling('Circ2\', 'Q1\', '08-26-19\', 0, 1113);
% Analyze_QP_Tunneling('Circ2\', 'Q3\', '08-26-19\', 0, 2536);
% Analyze_QP_Tunneling('Circ2\', 'Q3\', '08-26-19\', 2537, 2619);
% Analyze_QP_Tunneling('Circ2\', 'Q3\', '09-06-19\', 0, 499);
% Analyze_QP_Tunneling('Circ2\', 'Q3\', '09-06-19\', 500, 793);
% Analyze_QP_Tunneling('Circ2\', 'Q4\', '08-26-19\', 0, 99);
% Analyze_QP_Tunneling('Circ1\', 'Q2\', '08-30-19\', 0, 1242);
Analyze_QP_Tunneling('Circ2\', 'Q4\', '09-06-19\', 0, 734);
Analyze_QP_Tunneling('Circ2\', 'Q4\', '09-07-19\', 0, 264);

% plotFit(1/0.28e-3, 2e-5);
% plotFit(1/13e-3, 9e-4);
% plotFit(1/1e-3, 5e-5);
% plotFit(1/1.5e-3, 4e-4);

function [] = plotFit(gamma, A0)
    f = 10.^[-1:0.01:3.5]; 
    A = gamma/4*A0;
    theoryplot = gamma*A./(4*pi^2*f.^2+gamma^2/4);
    figure(111); hold on; 
    plot(f,theoryplot, 'DisplayName', ['fit ',num2str(1e3/gamma),'ms'])
end