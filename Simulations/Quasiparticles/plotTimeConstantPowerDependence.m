function plotTimeConstantPowerDependence
%plotTimeConstantPowerDependence Plot the time constant bias dependence.

r = 1e-10; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 1e-3; % trapping rate in units of 1 / \tau_0
V = 1.5:.5:20; % in units of \Delta
Tph = 0.050; % K
tspan = [-5000, 5000]; % in units of \tau_0

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
P = NaN(size(V));
for k = 1:length(V)
    % [t, ~, ~, ~, n_qp, ~, ~, P(k)]  = noTrapping0DModel(r, V(k), Tph, tspan);
    [t, ~, ~, ~, n_qp, ~, ~, P(k)]  = simpleTrapping0DModel(r, c, V(k), Tph, tspan);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = extractTimeConstants(t, n_qp, false);
end

figure
errorbar(V, tau_p, err_p, '.', 'Color', [0.8500    0.3250    0.0980],...
    'LineWidth', 3)
hold on
errorbar(V, tau_r, err_r, '.', 'Color', [0    0.4470    0.7410],...
    'LineWidth', 3)
xlabel('Bias Voltage (\Delta/e)', 'FontSize', 14)
ylabel('Time Constant (\tau_0)', 'FontSize', 14)
legend('poisoning', 'recovery')
title({'Time Constants',...
      ['injection rate r = ', num2str(r), '/\tau_0, ',...
       'trapping rate c = ', num2str(c), '/\tau_0']})
set(gca, 'yscale', 'Log')
grid on
grid minor
axis tight

figure
errorbar(P, tau_p, err_p, '.', 'Color', [0.8500    0.3250    0.0980],...
    'LineWidth', 3)
hold on
errorbar(P, tau_r, err_r, '.', 'Color', [0    0.4470    0.7410],...
    'LineWidth', 3)
xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
ylabel('Time Constant (\tau_0)', 'FontSize', 14)
legend('poisoning', 'recovery')
title({'Time Constants',...
      ['injection rate r = ', num2str(r), '/\tau_0, ',...
       'trapping rate c = ', num2str(c), '/\tau_0']})
set(gca, 'yscale', 'Log')
grid on
grid minor
axis tight

end