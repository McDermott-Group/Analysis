function plotPowerDependences
%plotPowerDependences Plot the stationary (equilibrium)
%quasiparticle concentration vs injection power/voltage, as well
%as the time that it takes to achieve the equilibrium vs injection power/
%voltage.

Tph = 0.050; % K
V = 1:1:20; % in units of \Delta
r = 0.01; % in units of 1 / \tau_0 (assuming n_{qp} in units of n_{cp})
tspan = [-10, 0]; % in units of tau0

P = nan(size(V));
n_qp_eq = nan(size(V));
t_qp_eq = nan(size(V));
for k = 1:length(V)
    [~, ~, ~, ~, n_qp, ~, P(k)] = noTrapping0DModel(r, V(k), Tph, tspan);
    n_qp_eq(k) = n_qp(end);
    [~, idx] = max(n_qp);
    t_qp_eq(k) = n_qp(idx);
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
xlabel('Injection Voltage (eV/\Delta)', 'FontSize', 14)
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
ylabel('Poisoning Equilibriation Time (\tau_0)', 'FontSize', 14)
title({'Poisoning Equilibriation Time',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

figure
plot(V, t_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Voltage (eV/\Delta)', 'FontSize', 14)
ylabel('Poisoning Equilibriation Time (\tau_0)', 'FontSize', 14)
title({'Poisoning Equilibriation Time',...
    ['Injection Rate Constant r = ', num2str(r), '/\tau_0']})
grid on
grid minor
axis tight

end