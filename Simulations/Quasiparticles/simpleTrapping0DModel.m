function [t, e, n, f, n_qp, n_qp_T, r_qp, P] = ...
    simpleTrapping0DModel(r, c, V, Tph, tspan)
%simpleTrapping0DModel Simplest model with quasiparticle traps.
% [t, e, n, f, n_qp, n_qp_T, r_qp, P] = noTrapping0DModel(r, V, Tph, tspan)
%   computes the quasiparticle dynamics in a simple model without any
%   quasiparticle traps.
%   The input parameters:
%      r is the injection rate in 1/\tau_0, assuming n and n_qp in units
%      of n_{cp},
%      c is the trapping (capture) rate in 1/\tau_0, assuming n and n_qp
%      in units of n_{cp},
%      V is the applied voltage in \Delta,
%      Tph is phonon temperature in K,
%      tspan should be in form of [ti, tf] in units of \tau_0 where ti is
%      the inital time and tf is the final time of the integration (it is
%      assumed quasiparticles are injected at t < 0 and relaxation happens
%      at t > 0).
%   The output parameters:
%      t is time in \tau_0,
%      e are energies of the bins in \Delta,
%      n are the quasiparticle densities in n_{cp}, size(n) == [length(t),
%      length(e)],
%      f are the occupational numbers, size(f) == size(n),
%      n_qp is the non-equilibrium quasiparticle density,
%      length(n_qp) == length(t),
%      n_qp_T is the thermal quasiparticle density
%      length(n_qp_T) == lenght(t),
%      r_qp is the total injection rate in 1 / \tau_0, assuming n and n_qp
%      in units n_{cp}, single number,
%      P is the total injected power in \Delta / \tau_0, assuming n and
%      n_qp in units n_{cp}, single number.
% Constants.
kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

% Tc = 1.2; % K (aluminum critical temperature)
delta = 0.18e-3; % eV (aluminum superconducting gap)

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta

% Tc = kB * Tc / delta;
Tc = 1 / 1.764; % \delta/(K_B * T_c) = 1.764 at T = 0 K BCS result

% Number of energy bins.
N = 500;
% Maximum energy.
max_e = max(V) + 2;

% Assign energies to the bins. Non-uniform energy separation is
% implemented. For a uniform distribution set alpha to a small positive
% value.
alpha = 2.5;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

[ej, ei] = meshgrid(e);

% Matrix elements defined by Eq. (6) and Eq. (7) in J. M. Martinis et al.,
% Phys. Rev. Lett. 103, 097002 (2009), as well as Eqs. (C3) - (C6) in 
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
Gs = (ei - ej).^2 .* (1 - 1 ./ (ei .* ej)) .* rho(ej) ./...
    abs(exp(-(ei - ej) / Tph) - 1) / Tc^3 .*...
    (ones(length(e), 1) * de');
Gs(~isfinite(Gs)) = 0;

Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) ./...
    abs(exp(-(ei + ej) / Tph) - 1) / Tc^3;

% Shorthand notation.
rho_de =  rho(e) .* de;

% Initial condition n0. n is the deviation from the thermal equlibrium for
% a given temperature, i.e. the solution is the non-quilibrium correction.
% The total number of quasiparticle is actually
% n^{injected}_{qp} + n^T_{qp}(T), or n_qp + n_qp_T using the notation of
% this function.
n0 = zeros(size(rho_de));

% It is assmumed that the injection is happening at t < 0 and
% the relaxation at t > 0.
% Injection voltage.
R = zeros(size(e));
if length(V) == 2
    % Injection into an energy band.
    R(V(1) < e & e <= V(2)) = r;
else
    % Realistic injection.
    R(e <= V) = r;
end
R = R .* rho_de;

Gt = c * Gtrapping(e, Tph, Tc);

% Solve the ODE.
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n, rho_de, Gs, Gr, R, Gt),...
    tspan, n0, options);

% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipatical density.
n_qp = 2 * sum(n, 2);

% Thermal quasipatical density.
n_qp_T = 2 * sum(rho_de ./ (exp(e / Tph) + 1));

% Bins indices the quasiparticles are injected to.
if length(V) == 2
    indices = V(1) < e & e < V(2);
else
    indices = e < V;
end

% Total injection rate.
r_qp = r * sum(rho_de(indices));

% Injection power.
P = r * sum(e(indices) .* rho_de(indices));

end

function normalized_density = rho(e)
    normalized_density = e ./ sqrt(e.^2 - 1);
end

function Gt = Gtrapping(e, Tph, Tc)
    de = 1e-5; % in units of \Delta
    e_subgap = de/2:de:1-de/2;
    [e_subgap, ei] = meshgrid(e_subgap, e);

    Gt = (ei - e_subgap).^2 ./ abs(exp(-(ei - e_subgap) / Tph) - 1) / Tc^3;
    Gt(~isfinite(Gt)) = 0;
    Gt = sum(Gt, 2) * de;
end

function ndot = quasiparticleODE(t, n, rho_de, Gs, Gr, R, Gt)
    Gs = Gs .* (1 - ones(length(n), 1) * (n ./ rho_de)');
    if t > 0
        R = 0;
    end
    % Eq. (2) in J. Wenner and J. M. Martinis: Erratum to
    % M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
    % http://web.physics.ucsb.edu/~martinisgroup/papers/Lenander2010erratum.pdf
    % The last term is due to quasiparticle trapping assuming a simple
    % model such that a quasiparticle is trapped if it scatters below
    % the gap.
    ndot = Gs' * n - sum(Gs, 2) .* n - 2 * n .* (Gr * n) + R - Gt .* n;
end