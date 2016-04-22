function fit2NIS04212016
%FIT2NIS04212016 Fitting to FinalNIS04212016forSimulations data set.

data = load('NIS04212016.mat');

% V_p_f = data.FarTrapPoisoning(:, 1);
E_p_f = data.FarTrapPoisoning(:, 2);
tau_p_f = data.FarTrapPoisoning(:, 3);
nqp_p_f = data.FarTrapPoisoning(:, 4);

% V_r_f = data.FarTrapRecovery(:, 1);
E_r_f = data.FarTrapRecovery(:, 2);
tau_r_f = data.FarTrapRecovery(:, 3);
nqp_r_f = data.FarTrapRecovery(:, 4);

% V_p_n = data.NearTrapPoisoning(:, 1);
E_p_n = data.NearTrapPoisoning(:, 2);
tau_p_n = data.NearTrapPoisoning(:, 3);
nqp_p_n = data.NearTrapPoisoning(:, 4);

% V_r_n = data.NearTrapRecovery(:, 1);
E_r_n = data.NearTrapRecovery(:, 2);
tau_r_n = data.NearTrapRecovery(:, 3);
nqp_r_n = data.NearTrapRecovery(:, 4);

r = 1.5e-7; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = .008; % trapping rate in units of 1 / \tau_0
V = [1.5, 2, 3, 4, 5, 10, 20]; % in units of \Delta
Tph = 0.051; % K
tspan = [-500, 500]; % in units of \tau_0

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
nqp = NaN(size(V));
P = NaN(size(V));
for k = 1:length(V)
    [t, ~, ~, ~, n_qp, ~, ~, P(k)] = simpleTrapping0DModel(Tph, tspan, V(k), r, c);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = extractTimeConstants(t, n_qp, false);
    nqp(k) = max(n_qp);
end

nqp = 4e6  * nqp; % n_{cp} for aluminum is (4e-6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)

tau0 = 5 * .438; % us, \tau_0 for aluminum from S. B. Kaplan et al.,
              % Phys. Rev. B 14, 4854 (1976)
tau_p = tau0 * tau_p;
tau_r = tau0 * tau_r;

scrsz = get(0, 'ScreenSize');
figure('Position', [.1 .1 1.5 .8] * scrsz(4));
subplot(1, 2, 1)
hold on
plot(E_p_f, tau_p_f, 'o',...
     E_r_f, tau_r_f, 's',...
     E_p_n, tau_p_n, '^',...
     E_r_n, tau_r_n, 'v', 'MarkerSize', 10, 'LineWidth', 2)
errorbar(V, tau_p, err_p, '.', 'LineWidth', 2)
errorbar(V, tau_r, err_r, '.', 'LineWidth', 2)
set(gca, 'xscale', 'Log')
hold off
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Time Constant (\mu s)', 'FontSize', 14)
% legend({'poisoning, far trap', 'recovery, far trap',...
%     'poisoning, near trap', 'recovery, near trap'})
title('Time Constants (Experiment vs Simulation)')
axis tight
grid on
 
subplot(1, 2, 2)
hold on
loglog(E_p_f, nqp_p_f, 'o',...
       E_r_f, nqp_r_f, 's',...
       E_p_n, nqp_p_n, '^',...
       E_r_n, nqp_r_n, 'v',...
       V, nqp, '*',...
       'MarkerSize', 10, 'LineWidth', 2)
hold off
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'poisoning, far trap', 'recovery, far trap',...
    'poisoning, near trap', 'recovery, near trap'},...
    'Location', 'NorthWest')
title('Quasiparticle Steady-State Density (Experiment vs Simulation)')
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
xlim([1 100])
axis tight
grid on
end