function steadyStateModelPowerLowExponent
%steadyStateModelPowerLowExponent Explore the quasiparticle steady-state
% density at different injection rates (resistances).

r_direct = logspace(-2, 2, 5) * 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = .5; % dimensionless
c = [0, 0.005, .01, .02, .05, .1, .15, .2]; % dimensionless
vol = 2.6e+04; % um^3
V = linspace(1.001, 1.5, 20); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nqp = NaN(length(V), length(r_direct));
P = NaN(size(nqp));
power = NaN(size(r_direct));
power_mean = NaN(size(c));
power_err = NaN(size(c));
for kc = 1:length(c)
    for krqp = 1:length(r_direct)
        rqpsim = r_direct(krqp);
        parfor kV = 1:length(V)
            Vsim = V(kV);
            csim = c(kc);
            [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
                twoRegionSteadyStateModelOptimized(Tph, tspan,...
                Vsim, rqpsim, r_phonon, csim, vol, N);
            nqp(kV, krqp) = max(n_qp);
            fprintf('*')
        end
        p = polyfit(log(P(:, krqp)), log(nqp(:, krqp)), 1);
        power(krqp) = p(1);
    end
    power_mean(kc) = mean(power);
    power_err(kc) = std(power) / sqrt(length(power));
end

h = figure;
errorbar(c, power_mean, power_err, '.', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('Trapping Strength c_{tr}', 'FontSize', 14)
ylabel('Power-Low Exponent \alpha', 'FontSize', 14)
title('n_{qp}^{eq} \propto P^{\alpha}')
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimPowerLowExponent.pdf')
end