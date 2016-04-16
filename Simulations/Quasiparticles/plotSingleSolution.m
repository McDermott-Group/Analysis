function plotSingleSolution
%PLOTSINGLESOLUTION Plot the quasiparticle dynamics solution. 
%   Detailed explanation goes here

Tph = 0.050; % K
V = 1.2; % in units of \Delta
r = 0.1;
tspan = [-10, -9]; % in units of tau0
[t, e, ~, f, n_qp] = simplest0DModel(r, V, Tph, tspan);

figure
plot(t, n_qp, 'LineWidth', 3)
hold on
xlabel('Time (t/\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
title({'Quasipaticle Dynamics',...
       '(injection at t < 0, relaxation at t > 0)'})
grid on
grid minor
axis tight

figure
plotSmooth(t, e, f)
xlabel('Time (t/\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, relaxation at t > 0)'})

end