function plotSingleSolution
%PLOTSINGLESOLUTION Generate the quasiparticle dynamics plots.

r = 1e-8; % in units of 1 / \tau_0
% r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1));
                                      %(assuming n_{qp} in units of n_{cp})
c = 1e-3; % trapping rate in units of 1 / \tau_0
V = 1.1; %[2.8, 3]; % in units of \Delta
Tph = .05; % K
tspan = [-5000, 5000]; % in units of \tau_0

% [t, e, ~, f, n_qp] = noTrapping0DModel(r, V, Tph, tspan);
[t, e, ~, f, n_qp] = simpleTrapping0DModel(r, c, V, Tph, tspan);

figure
plot(t, n_qp, 'LineWidth', 3)
hold on
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
title({'Quasipaticle Dynamics',...
       '(injection at t < 0, recovery at t > 0)'})
grid on
grid minor
axis tight

figure
plotSmooth(t, e, f)
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, recovery at t > 0)'})

extractTimeConstants(t, n_qp, true);

end