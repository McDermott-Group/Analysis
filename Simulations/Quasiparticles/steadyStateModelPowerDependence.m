function steadyStateModelPowerDependence
%steadyStateModelPowerDependence Exploring power dependence at different
%biases.

r_direct = [.01, .1, 1, 10, 100] * 1e-07; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = [1, .1, .01, .001, .0001]; % dimensionless
c = 0*2.5e-02; % dimensionless
vol = 5.000e+03; % um^3
V = 1.1:1:30; % in units of \delta

Tph = 0.051; % K
tspan = [-300, 0]; % in units of \tau_0

% Number of the energy bins.
N = 500;

nqp = NaN(length(V), length(r_direct));
P = NaN(size(nqp));
for krqp = 1:length(r_direct)
    rqpsim = r_direct(krqp);
    rphsim = r_phonon(krqp);
    parfor kV = 1:length(V)
        Vsim = V(kV);
        [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
            twoRegionSteadyStateModel(Tph, tspan,...
            Vsim, rqpsim, rphsim, c, vol, N, false);
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
    legends{k} = ['r_{qp} = ', num2str(r_direct(k), '%.2e'), '\Delta/e'];
end
legend(legends, 'Location', 'SouthEast')
title({'Power Dependence at Different Bias Voltages',...
    ['r_{ph} = ', num2str(r_phonon, '%.2e'), ', ',...
     'c_{trap} = ', num2str(c, '%.2e'), ', '...
     'volume = ', num2str(vol, '%.2e'), '\mu{m}^3']})
axis tight
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
grid on

end