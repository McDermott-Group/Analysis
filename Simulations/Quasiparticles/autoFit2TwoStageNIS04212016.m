function autoFit2TwoStageNIS04212016
%autoFit2TwoStageTrapNoTrap Fitting to TrapNoTrap dataset.

% r_direct = 1.676e-05; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
% r_phonon = 4.876e-02; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
% c = 1.603e-02; % trapping rate in units of 1 / \tau_0
r_direct = 8.323e-06; r_phonon = 5.018e-03; c = 1.639e-02;

Tph = 0.051; % K
tspan = [-310, -10]; % in units of \tau_0

% Number of the energy bins.
N = 100;

% Time domain data fit.
data = load('NIS04212016.mat');
E_p_n = data.NearTrapPoisoning(:, 2);
nqp_p_n = data.NearTrapPoisoning(:, 4);

E_r_n = data.NearTrapRecovery(:, 2);
nqp_r_n = data.NearTrapRecovery(:, 4);

V = [E_p_n; E_r_n];
nqp = [nqp_p_n; nqp_r_n];

options = optimset('Display', 'iter', 'MaxIter', floor(2000 / N), 'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, nqp, N),...
    [r_direct, r_phonon, c], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), ';'])
end

function error = simulations(x, Tph, tspan, V, nqp, N)
    indices = (V > 1) & (nqp > 0);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    for k = 1:length(V)
        [~, ~, ~, ~, n_qp] = twoStageQuasi0DModel(Tph, tspan,...
            V(k), x(1), x(2), x(3), 1e5, N, false);
        nqp_sim(k) = n_qp(end);
    end
    error = sum(log(nqp_sim ./ nqp).^2);
end