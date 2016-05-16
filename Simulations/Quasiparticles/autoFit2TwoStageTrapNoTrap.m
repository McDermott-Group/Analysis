function autoFit2TwoStageTrapNoTrap
%autoFit2TwoStageTrapNoTrap Fitting to TrapNoTrap dataset.

% No traps.
r_direct = 1.676e-05; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 4.876e-02; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 1.603e-02; % trapping rate in units of 1 / \tau_0
vol = 3.044e+04; % um^3

% With traps.
% r_direct = 4.875e-06; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
% r_phonon = 1.138e-02; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
% c = 1.591e-02; % trapping rate in units of 1 / \tau_0
% vol = 2.585e+05; % um^3

Tph = 0.051; % K
tspan = [-310, -10]; % in units of \tau_0

% Number of the energy bins.
N = 50;

delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('TrapNoTrap.mat');

No traps.
V = data.NoTrap(:, 5) / delta;
P = data.NoTrap(:, 6);
nqp = data.NoTrap(:, 8) - min(data.NoTrap(:, 8));

% With traps.
% V = data.Trap(:, 5) / delta;
% P = data.Trap(:, 6);
% nqp = data.Trap(:, 8) - min(data.Trap(:, 8));

options = optimset('Display', 'iter', 'MaxIter', floor(2000 / N), 'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, P, nqp, N),...
    [r_direct, r_phonon, c, vol], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(x(4), '%.3e'), '; '])
end

function error = simulations(x, Tph, tspan, V, P, nqp, N)
    indices = (V > 1) & (nqp > 0);
    P = P(indices);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    P_sim = NaN(size(V));
    for k = 1:length(V)
        [~, ~, ~, ~, n_qp, ~, P_sim(k)] = twoStageQuasi0DModel(Tph, tspan,...
            V(k), x(1), x(2), x(3), x(4), N, false);
        nqp_sim(k) = n_qp(end);
    end
    error = sum(log(nqp_sim ./ nqp).^2 + log(P_sim ./ P).^2);
end