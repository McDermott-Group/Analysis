function plotTimeConstantBiasDependence
%PLOTSINGLESOLUTION Plot the time constant bias dependence.

Tph = 0.050; % K
V = 3:1:10; % in units of \Delta
r = 0.0001;

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
for k = 1:length(V)
    [t, ~, ~, ~, n_qp] = simplest0DModel(r, V(k), Tph);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = extractTimeConstants(t, n_qp);
end

figure
errobar(V, tau_p, err_p, '-o', 'LineWidth', 3)
hold on
errobar(V, tau_r, err_r, '-o', 'LineWidth', 3)
xlabel('Bias Voltage (eV/\Delta)', 'FontSize', 14)
ylabel('Time Constant (\tau/\tau_0)', 'FontSize', 14)
title('Time Constants')
legend('poisonning', 'relaxation')
grid on
grid minor
axis tight

end