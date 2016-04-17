function plotSingleSolution
%PLOTSINGLESOLUTION Generate the quasiparticle dynamics plots.

Tph = .1; % K
V = [2.8, 3]; % in units of \Delta
r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1));

tspan = [-30000, 0]; % in units of tau0

[t, e, ~, f, n_qp] = noTrapping0DModel(r, V, Tph, tspan);

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