function plotSingleSolution
%PLOTSINGLESOLUTION Generate the quasiparticle dynamics plots.

r = 1e-9; % 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1));
V = 2; %[2.8, 3]; % in units of \Delta
Tph = .05; % K
tspan = [-30000, 50000]; % in units of \tau_0

[t, e, ~, f, n_qp] = noTrapping0DModel(r, V, Tph, tspan);

figure
plot(t, n_qp, 'LineWidth', 3)
hold on
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
title({'Quasipaticle Dynamics',...
       '(injection at t < 0, relaxation at t > 0)'})
grid on
grid minor
axis tight

figure
plotSmooth(t, e, f)
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, relaxation at t > 0)'})

end