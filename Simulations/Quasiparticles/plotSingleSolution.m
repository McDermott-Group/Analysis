function plotSingleSolution
%PLOTSINGLESOLUTION Generate the quasiparticle dynamics plots.

% r = 1e-8; % in units of 1 / \tau_0
% r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1));
%                                       %(assuming n_{qp} in units of n_{cp})
% c = 1; % trapping rate in units of 1 / \tau_0
% V = 2; %[2.8, 3]; % in units of \Delta
% Tph = .050; % K
% tspan = [-10000, 10000]; % in units of \tau_0

r_direct = 8.323e-06; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 5.018e-03; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 1.639e-02; % trapping rate in units of 1 / \tau_0
vol = 5e+04; % um^3

N = 125;

Tph = 0.051; % K
tspan = [-200, 200]; % in units of \tau_0

V = 2.5;

% [t, e, n, f, n_qp] = twoStageQuasi0DModel(Tph, tspan, V, r_direct, r_phonon, c, vol, N, true);
clear twoStageTimeDomainQuasi0DModel
[t, e, ~, f, n_qp] = twoStageTimeDomainQuasi0DModel(Tph, tspan, V, r_direct, r_phonon, c, vol, N);

figure
plot(t, n_qp, 'LineWidth', 3)
hold on
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp} (\mu m^{-3})', 'FontSize', 14)
title({'Quasipaticle Dynamics',...
       '(injection at t < 0, recovery at t > 0)'})
grid on
grid minor
axis tight

figure
plotSmooth(t, e, f)
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, recovery at t > 0)'})

% figure
% semilogy(e, n(end, :), e, f(end, :), 'LineWidth', 3)
% xlabel('Energy (\Delta)', 'FontSize', 14)
% ylabel('n(\epsilon), f(\epsilon)', 'FontSize', 14)
% legend('n(\epsilon)', 'f(\epsilon)')
% axis tight
% grid on

extractTimeConstants(t, n_qp, true);

end