function autoFit2NIS04212016
%autoFit2NIS04212016 Fitting to the NIS04212016 dataset.

r_direct = 8.263e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 2.223e-01; % dimensionless
c = 9.414e-02; % dimensionless
vol = 5e3; % um^3

Tph = 0.051; % K
tspan = [-310, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Time domain data fit.
data = load('NIS04212016.mat');
E_p_n = data.NearTrapPoisoning(:, 2);
nqp_p_n = data.NearTrapPoisoning(:, 4);

E_r_n = data.NearTrapRecovery(:, 2);
nqp_r_n = data.NearTrapRecovery(:, 4);

V = [E_p_n; E_r_n];
nqp = [nqp_p_n; nqp_r_n];

options = optimset('Display', 'iter', 'MaxIter', floor(10),...
    'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, nqp, N),...
    [r_direct, r_phonon, c, vol], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(x(4), '%.3e'), ';'])
end

function error = simulations(x, Tph, tspan, V, nqp, N)
    indices = (V > 1) & (nqp > 0);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    r_qp = x(1);
    r_ph = x(2);
    c = x(3);
    vol = x(4);
    parfor k = 1:length(V)
        [~, ~, ~, ~, n_qp] = twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V(k), r_qp, r_ph, c, vol, N);
        nqp_sim(k) = n_qp(end);
    end
    error = sum(log(nqp_sim ./ nqp).^2);
end