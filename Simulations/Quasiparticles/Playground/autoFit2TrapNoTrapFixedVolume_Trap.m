function autoFit2TrapNoTrapFixedVolume_Trap
%autoFit2TrapNoTrapFixedVolume_Trap Fitting to the TrapNoTrap dataset.

With traps.
r_direct = 1.649e-04; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 8.646e-05; % dimensionless
c = 1.021e-02; % dimensionless
vol = 5e3; % um^3

Tph = 0.051; % K
tspan = [-310, -10]; % in units of \tau_0

% Number of the energy bins.
N = 150;

delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('TrapNoTrap.mat');

% With traps.
V = data.Trap(:, 5) / delta;
P = data.Trap(:, 6);
nqp = data.Trap(:, 8) - min(data.Trap(:, 8));

options = optimset('Display', 'iter', 'MaxIter', floor(5000 / N), 'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, P, nqp, vol, N),...
    [r_direct, r_phonon, c], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(vol, '%.3e'), '; '])
end

function error = simulations(x, Tph, tspan, V, P, nqp, vol, N)
    indices = (V > 1) & (nqp > 0);
    P = P(indices);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    P_sim = NaN(size(V));
    for k = 1:length(V)
        [~, ~, ~, ~, n_qp, ~, P_sim(k)] = ...
            phononMediatedPoisoningEquilibriumModel(Tph, tspan,...
            V(k), x(1), x(2), x(3), vol, N, false);
        nqp_sim(k) = n_qp(end);
    end
    error = sum(log(nqp_sim ./ nqp).^2 + log(P_sim ./ P).^2);
end