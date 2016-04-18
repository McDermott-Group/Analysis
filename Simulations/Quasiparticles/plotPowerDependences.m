function plotPowerDependences
%plotPowerDependences Plot the stationary (equilibrium)
%quasiparticle concentration vs injection power/voltage, as well
%as the time that it takes to achieve the equilibrium vs injection power/
%voltage.

r = 1e-9; % in units of 1 / \tau_0 (assuming n_{qp} in units of n_{cp})
V = 1:.25:20; % in units of \Delta
Tph = 0.050; % K
tspan = [-30000, 0]; % in units of \tau_0

P = nan(size(V));
n_qp_eq = nan(size(V));
t_qp_eq = nan(size(V));
for k = 1:length(V)
    [t_qp, ~, ~, ~, n_qp, ~, ~, P(k)] = noTrapping0DModel(r, V(k), Tph, tspan);
    n_qp_eq(k) = n_qp(end);
    [~, idx] = max(n_qp);
    t_qp_eq(k) = t_qp(idx) - t_qp(1);
end

figure
plot(P, n_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp}^{eq}/n_{\rm cp}', 'FontSize', 14)
title({'Stationary (Equilibrium) Quasiparticle Concentration',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

figure
plot(V, n_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Voltage (\Delta/e)', 'FontSize', 14)
ylabel('n_{\rm qp}^{eq}/n_{\rm cp}', 'FontSize', 14)
title({'Stationary (Equilibrium) Quasiparticle Concentration',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

figure
plot(P, t_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
ylabel('Poisoning Equilibration Time (\tau_0)', 'FontSize', 14)
title({'Poisoning Equilibration Time',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

figure
plot(V, t_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Voltage (\Delta/e)', 'FontSize', 14)
ylabel('Poisoning Equilibration Time (\tau_0)', 'FontSize', 14)
title({'Poisoning Equilibration Time',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

end