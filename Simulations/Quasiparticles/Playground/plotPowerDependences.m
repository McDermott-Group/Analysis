function plotPowerDependences
%plotPowerDependences Plot the stationary (equilibrium)
%quasiparticle concentration vs injection power/voltage, as well
%as the time that it takes to achieve the equilibrium vs injection power/
%voltage.

r = 1e-8; % in units of 1 / \tau_0 (assuming n_{qp} in units of n_{cp})
c = 1e-3; % trapping rate in units of 1 / \tau_0
V = [1, 1.01, 1.03, 1.07, 1.1, 1.25, 1.25:.5:10]; % in units of \Delta
Tph = 0.050; % K
tspan = [-20000, 0]; % in units of \tau_0

P_no = nan(size(V));
n_qp_eq_no = nan(size(V));
t_qp_eq_no = nan(size(V));
P_tr = nan(size(V));
n_qp_eq_tr = nan(size(V));
t_qp_eq_tr = nan(size(V));
for k = 1:length(V)
    [t_qp, ~, ~, ~, n_qp, ~, ~, P_no(k)] = noTrapping0DModel(r, V(k), Tph, tspan);
    n_qp_eq_no(k) = n_qp(end);
    [~, idx] = max(n_qp);
    t_qp_eq_no(k) = t_qp(idx) - t_qp(1);
    [t_qp, ~, ~, ~, n_qp, ~, ~, P_tr(k)] = simpleTrapping0DModel(r, c, V(k), Tph, tspan);
    n_qp_eq_tr(k) = n_qp(end);
    [~, idx] = max(n_qp);
    t_qp_eq_tr(k) = t_qp(idx) - t_qp(1);
end

figure
plot(P_no, n_qp_eq_no, P_tr, n_qp_eq_tr, 'LineWidth', 3)
hold on
xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp}^{eq}/n_{\rm cp}', 'FontSize', 14)
title({'Stationary (Equilibrium) Quasiparticle Concentration',...
    ['injection rate r = ', num2str(r), '/\tau_0, ',...
     'trapping rate c = ', num2str(c), '/\tau_0']})
legend('no trapping', 'with trapping')
grid on
grid minor
axis tight

figure
plot(V, n_qp_eq_no, V, n_qp_eq_tr, 'LineWidth', 3)
hold on
xlabel('Injection Voltage (\Delta/e)', 'FontSize', 14)
ylabel('n_{\rm qp}^{eq}/n_{\rm cp}', 'FontSize', 14)
title({'Stationary (Equilibrium) Quasiparticle Concentration',...
    ['injection rate r = ', num2str(r), '/\tau_0, ',...
     'trapping rate c = ', num2str(c), '/\tau_0']})
legend('no trapping', 'with trapping')
grid on
grid minor
axis tight

% figure
% plot(P_no, t_qp_eq_no, P_tr, t_qp_eq_tr, 'LineWidth', 3)
% hold on
% xlabel('Injection Power (n_{\rm cp}\Delta/\tau_0)', 'FontSize', 14)
% ylabel('Poisoning Equilibration Time (\tau_0)', 'FontSize', 14)
% title({'Stationary (Equilibrium) Quasiparticle Concentration',...
%     ['injection rate r = ', num2str(r), '/\tau_0, ',...
%      'trapping rate c = ', num2str(c), '/\tau_0']})
% grid on
% grid minor
% axis tight
% 
% figure
% plot(V, t_qp_eq_no, V, t_qp_eq_tr, 'LineWidth', 3)
% hold on
% xlabel('Injection Voltage (\Delta/e)', 'FontSize', 14)
% ylabel('Poisoning Equilibration Time (\tau_0)', 'FontSize', 14)
% title({'Stationary (Equilibrium) Quasiparticle Concentration',...
%     ['injection rate r = ', num2str(r), '/\tau_0, ',...
%      'trapping rate c = ', num2str(c), '/\tau_0']})
% legend('no trapping', 'with trapping')
% grid on
% grid minor
% axis tight

end