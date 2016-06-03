function timeDomainModelFit2NIS04212016
%timeDomainModelFit2NIS04212016 Fitting to the FinalNIS04212016forSimulations
% dataset using the two-point time domain quasi-0D model.

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

r_direct = 8.263e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 2.223e-01; % dimensionless
c = 9.414e-02; % dimensionless
vol = 5e3; % um^3

V = [E_p_n; E_r_n; 3.5; 4.2; 6.7]; % in units of \Delta
Tph = 0.051; % K
tspan = [-200, 200]; % in units of \tau_0

% Number of energy bins.
N = 150;

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
nqp = NaN(size(V));
for k = 1:length(V)
    clear twoRegionTimeDomainModel
    [t, ~, ~, ~, n_qp] = ...
        twoRegionTimeDomainModel(Tph, tspan, V(k),...
        r_direct, r_phonon, c, vol, N);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = ...
        extractTimeConstants(t, n_qp, false);
    nqp(k) = max(n_qp);
    fprintf('*')
end
fprintf('\n')

tau0 = .438; % us, \tau_0 for aluminum from S. B. Kaplan et al.,
              % Phys. Rev. B 14, 4854 (1976)
tau_p = tau0 * tau_p;
tau_r = tau0 * tau_r;

F = median(tau_p_n ./ tau_p(1:length(tau_p_n)));

tau_p = F * tau_p;
tau_r = F * tau_r;

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
title({'Time Constants', ['r_{qp} = ', num2str(r_direct, '%.3e'),...
    ', r_{ph} = ', num2str(r_phonon, '%.3e'),...
    ', c_{tr} = ', num2str(c, '%.3e'),...
    ', F = ', num2str(F, '%.3f')]})
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

saveas(gca, 'NIS24062016.pdf', 'pdf')

end