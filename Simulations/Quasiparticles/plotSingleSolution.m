function plotSingleSolution
%PLOTSINGLESOLUTION Generate the quasiparticle dynamics plots.

% r = 1e-8; % in units of 1 / \tau_0
% r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1));
%                                       %(assuming n_{qp} in units of n_{cp})
% c = 1; % trapping rate in units of 1 / \tau_0
% V = 2; %[2.8, 3]; % in units of \Delta
% Tph = .050; % K
% tspan = [-10000, 10000]; % in units of \tau_0

r_direct = .09e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 1; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 0; % trapping rate in units of 1 / \tau_0
vol = 1e5; % um^- 3

N = 1000;

Tph = 0.051; % K
tspan = [-500, -1]; % in units of \tau_0

V = 4;

[t, e, n, f, n_qp] = twoStageQuasi0DModel(Tph, tspan, V, r_direct, r_phonon, c, vol, N, true);

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

figure
semilogy(e, n(end, :), e, f(end, :), 'LineWidth', 3)
xlabel('Energy (\Delta)', 'FontSize', 14)
ylabel('n(\epsilon), f(\epsilon)', 'FontSize', 14)
legend('n(\epsilon)', 'f(\epsilon)')
axis tight
grid on

% extractTimeConstants(t, n_qp, true);

end