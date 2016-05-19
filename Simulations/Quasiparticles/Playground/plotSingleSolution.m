function plotSingleSolution
%plotSingleSolution Quasiparticle dynamics plots.

r_direct = 8.323e-06; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 5.018e-03; % dimensionless
c = 1.639e-02; % dimensionless
vol = 5e+04; % um^3

N = 125;

Tph = 0.051; % K
tspan = [-200, 200]; % in units of \tau_0

V = 2.5;

% [t, e, n, f, n_qp] = ...
    % phononMediatedPoisoningEquilibriumModel(Tph, tspan, V,...
    % r_direct, r_phonon, c, vol, N, true);
clear phononMediatedPoisoningTimeDomainModel
[t, e, n, f, n_qp] = ...
    phononMediatedPoisoningTimeDomainModel(Tph, tspan, V,...
    r_direct, r_phonon, c, vol, N);

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
[Ind1, Ind2] = ndgrid(t, e);
hndl = surf(Ind1, Ind2, f);
set(gca, 'View', [0 90])
set(hndl, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong');
axis tight
colormap(jet)
colorbar
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, recovery at t > 0)'})

figure
n(n < 0) = NaN;
f(f < 0) = NaN;
semilogy(e, n(end, :), e, f(end, :), 'LineWidth', 3)
xlabel('Energy (\Delta)', 'FontSize', 14)
ylabel('n(\epsilon), f(\epsilon)', 'FontSize', 14)
legend('n(\epsilon)', 'f(\epsilon)')
axis tight
grid on

extractTimeConstants(t, n_qp, true);

end