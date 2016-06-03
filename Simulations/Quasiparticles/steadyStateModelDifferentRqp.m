function steadyStateModelDifferentRqp
%steadyStateModelDifferentRqp Explore the quasiparticle steady-state
% density at different injection rates (resistances).

r_direct = [.1, 1, 10, 100] * 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 5e-2; % dimensionless
c = 0; % dimensionless
vol = 5.000e+03; % um^3
V = 1.1:1:30; % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nqp = NaN(length(V), length(r_direct));
P = NaN(size(nqp));
for krqp = 1:length(r_direct)
    rqpsim = r_direct(krqp);
    parfor kV = 1:length(V)
        Vsim = V(kV);
        [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            Vsim, rqpsim, r_phonon, c, vol, N);
        nqp(kV, krqp) = max(n_qp);
        fprintf('*')
    end
end
fprintf('\n')
figure
hold on
for k = 1:length(r_direct)
    plot(P(:, k), nqp(:, k), 'MarkerSize', 10, 'LineWidth', 2)
end
xlabel('Injection Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legends = cell(0);
for k = 1:length(r_direct)
    legends{k} = ['r_{qp} = ', num2str(r_direct(k), '%.2e'), '/\tau_0'];
end
legend(legends, 'Location', 'SouthEast')
title(['r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', '...
     'v = ', num2str(vol, '%.2e'), ' \mu{m}^3'])
axis tight
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
grid on

saveas(gca, 'rqp.pdf', 'pdf')
end