function fit2NIS04212016
%FIT2NIS04212016 Fitting to FinalNIS04212016forSimulations data set.

data = load('NIS04212016.mat');

% V_p_f = data.FarTrapPoisoning(:, 1);
% E_p_f = data.FarTrapPoisoning(:, 2);
% tau_p_f = data.FarTrapPoisoning(:, 3);
% nqp_p_f = data.FarTrapPoisoning(:, 4);

% V_r_f = data.FarTrapRecovery(:, 1);
% E_r_f = data.FarTrapRecovery(:, 2);
% tau_r_f = data.FarTrapRecovery(:, 3);
% nqp_r_f = data.FarTrapRecovery(:, 4);

% V_p_n = data.NearTrapPoisoning(:, 1);
E_p_n = data.NearTrapPoisoning(:, 2);
tau_p_n = data.NearTrapPoisoning(:, 3);
nqp_p_n = data.NearTrapPoisoning(:, 4);

% V_r_n = data.NearTrapRecovery(:, 1);
E_r_n = data.NearTrapRecovery(:, 2);
tau_r_n = data.NearTrapRecovery(:, 3);
nqp_r_n = data.NearTrapRecovery(:, 4);

r = .5e-9; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 0.01; % trapping rate in units of 1 / \tau_0
d = 1;
V = [E_p_n; E_r_n; 3.5; 4.2; 6.7]; % in units of \Delta
% V = V(V > 3);
Tph = 0.051; % K
tspan = [-100, 100]; % in units of \tau_0

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
nqp = NaN(size(V));
P = NaN(size(V));
for k = 1:length(V)
    % [t, ~, ~, ~, n_qp, ~, ~, P(k)] = simpleTrapping0DModel(Tph, tspan, V(k), r, c);
    [t, ~, ~, ~, n_qp, ~, ~, P(k)] = simpleTrappingQuasi0DModel(Tph, tspan, V(k), r, c, d);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = extractTimeConstants(t, n_qp, false);
    nqp(k) = max(n_qp);
end

nqp = 4e6 * nqp; % n_{cp} for aluminum is (4e-6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)

F = 5;
tau0 = F * .438; % us, \tau_0 for aluminum from S. B. Kaplan et al.,
              % Phys. Rev. B 14, 4854 (1976)
tau_p = tau0 * tau_p;
tau_r = tau0 * tau_r;

scrsz = get(0, 'ScreenSize');
figure('Position', [.1 .1 1.5 .8] * scrsz(4));
subplot(1, 2, 1)
hold on
plot(E_p_n, tau_p_n, '^',...
     E_r_n, tau_r_n, 'v', 'MarkerSize', 10, 'LineWidth', 2)
errorbar(V, tau_p, err_p, '*', 'MarkerSize', 10, 'LineWidth', 2)
errorbar(V, tau_r, err_r, '*', 'MarkerSize', 10, 'LineWidth', 2)
set(gca, 'xscale', 'Log')
hold off
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Time Constant (\mu s)', 'FontSize', 14)
legend({'poisoning, near trap', 'recovery, near trap',...
    'poisoning simulation', 'recovery simulation'},...
    'Location', 'SouthWest')
title({'Time Constants', ['r = ', num2str(r, '%.2e'),...
    ', c = ', num2str(c, '%.2e'), ', d = ', num2str(d, '%.2e'),...
    ', F = ', num2str(F, '%.2f')]})
axis tight
grid on
 
subplot(1, 2, 2)
hold on
loglog(E_p_n, nqp_p_n, '^',...
       E_r_n, nqp_r_n, 'v',...
       V, nqp, '*',...
       'MarkerSize', 10, 'LineWidth', 2)
hold off
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'poisoning, near trap', 'recovery, near trap', 'simulation'},...
    'Location', 'NorthWest')
title('Quasiparticle Steady-State Density')
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
xlim([1 100])
axis tight
grid on

end