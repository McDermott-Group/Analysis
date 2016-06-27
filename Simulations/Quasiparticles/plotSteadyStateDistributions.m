function plotSteadyStateDistributions
%plotSingleSolution Quasiparticle dynamics plots.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct =1e-05; r_phonon = 5e-01; c = 0; vol = 2.6e+04; % um^3

N = 250;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

V = 1.5;

[~, e, n, f, ~, ~, ~, ~, ~, ~, ~, f_nis, n_nis] = ...
    twoRegionSteadyStateModelOptimized(Tph, tspan, V,...
    r_direct, r_phonon, c, vol, N);

h = figure;
n(n < 0) = NaN;
f(f < 0) = NaN;
n_nis(n_nis < 0) = NaN;
f_nis(n_nis < 0) = NaN;
semilogy(e, n_nis(end, :), e, f_nis(end, :),...
         e, n(end, :), e, f(end, :), 'LineWidth', 3)
xlabel('Energy (\epsilon=E/\Delta)', 'FontSize', 14)
ylabel('n/n_{cp}, f', 'FontSize', 14)
legend('n_{NIS}(\epsilon)', 'f_{NIS}(\epsilon)',...
       'n_{res}(\epsilon)', 'f_{res}(\epsilon)')
axis tight
xlim([1, max(V)])
grid on

c = 0.02;
[~, e, n_tr, f_tr, ~, ~, ~, ~, ~, ~, ~, f_nis_tr, n_nis_tr] = ...
    twoRegionSteadyStateModelOptimized(Tph, tspan, V,...
    r_direct, r_phonon, c, vol, N);

h = figure;
n_tr(n_tr < 0) = NaN;
f_tr(f_tr < 0) = NaN;
n_nis_tr(n_nis_tr < 0) = NaN;
f_nis_tr(n_nis_tr < 0) = NaN;
semilogy(e, f_nis(end, :), e, f_nis_tr(end, :),...
         e, f(end, :), e, f_tr(end, :), 'LineWidth', 3)
xlabel('Energy (\epsilon=E/\Delta)', 'FontSize', 14)
ylabel('Occupation Numbers f', 'FontSize', 14)
legend('f_{NIS}(\epsilon) (c=0)', ['f_{NIS}(\epsilon) (c=', num2str(c),')'],...
       'f_{res}(\epsilon) (c=0)', ['f_{res}(\epsilon) (c=', num2str(c),')'])
axis tight
xlim([1, max(V)])
grid on

savePDF(h, 'SimDistributions.pdf')
end