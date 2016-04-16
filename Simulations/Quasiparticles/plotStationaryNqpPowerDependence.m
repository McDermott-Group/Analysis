function plotStationaryNqpPowerDependence
%plotStationaryNqpPowerDependence Plot stationary (equilibrium)
%quasiparticle concentration vs injection power.
%   Detailed explanation goes here

Tph = 0.050; % K
V = 1:2:20; % in units of \Delta
r = 0.0001; % in units of 1 / \tau_0 (or n_{cp} / tau_0)
tspan = [-10, 0]; % in units of tau0

P = nan(size(V));
n_qp_eq = nan(size(V));
for k = 1:length(V)
    [~, ~, ~, ~, n_qp, ~, P(k)] = simplest0DModel(r, V(k), Tph, tspan);
    n_qp_eq(k) = n_qp(end);
end

figure
plot(P, n_qp_eq, 'LineWidth', 3)
hold on
xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp}^{eq}/n_{\rm cp}', 'FontSize', 14)
title({'Stationary (Equilibrium) Quasiparticle Concentration',...
    ['Injection Rate Constant r = ', num2str(r), ' n_{\rm cp}/\tau_0']})
grid on
grid minor
axis tight
end