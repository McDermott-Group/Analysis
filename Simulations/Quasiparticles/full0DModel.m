function [t, e, n, f, n_qp, n_qp_T, r_qp, P] = ...
    full0DModel(Tph, tspan, V, r, c)
%full0DModel Full 0D mode.
% [t, e, n, f, n_qp, n_qp_T, r_qp, P] = full0DModel(Tph, tspan, V, r, c)
%   computes the quasiparticle dynamics in a simple model without any
%   quasiparticle traps.
%   The input parameters:
%      Tph is phonon temperature in K,
%      tspan should be in form of [ti, tf] in units of \tau_0 where ti is
%      the inital time and tf is the final time of the integration (it is
%      assumed quasiparticles are injected at t < 0 and relaxation happens
%      at t > 0),
%      V is the applied voltage in \Delta,
%      r is the injection rate in 1/\tau_0, assuming n and n_qp in units
%      of n_{cp},
%      c is the trapping (capture) rate in 1/\tau_0, assuming n and n_qp
%      in units of n_{cp}.
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

% Debye temperature.
% TD = 428; % K
% TD = kB * TD/ delta; % in units of \Delta

% Number of energy bins.
N = 512;
% Maximum energy.
max_e = max(V) + 2;

% Assign the quasiparicle energies to the bins. Non-uniform energy
% separation is implemented. For a close to uniform distribution set alpha
% to a small positive value.
alpha = 2.5;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

% Short-hand notation.
rho_de = rho(e) .* de;

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gb = Gbreaking(e, de, Tph, Tc);

Gt = Gtrapping(e, Tph, Tc, c);

R = Injection(e, rho_de, V, r);

% Initial condition n0.
f0 = 1 ./ (exp(e / Tph) + 1);
n0 = rho_de .* f0;

% Debugging plot.
% figure
% plot(e, sum(Gs, 2), e, sum(Gr, 2), e, Gt, e, Gb, e, R, 'LineWidth', 2)
% xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
% ylabel('G (1/\tau_0)', 'FontSize', 14)
% legend({'G_{scattering}', 'G_{recombination}', 'G_{trapping}',...
%     'G_{pair-breaking}', 'R'})
% set(gca, 'yscale', 'Log')
% axis tight
% grid on

% Solve the ODE.
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gb, Gt, R), tspan, n0, options);

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

function abs_phonon_density = Np(e, Tph)
    abs_phonon_density = 1 ./ abs(exp(-e / Tph) - 1);
end

function [Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc)
% Scattering matrix defined by Eq. (6) in J. M. Martinis et al.,
% Phys. Rev. Lett. 103, 097002 (2009), as well as Eqs. (C3) and (C4) in 
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
% These equations come from Eq. (8) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).
    [ej, ei] = meshgrid(e);

    Gs = (ei - ej).^2 .* (1 - 1 ./ (ei .* ej)) .* rho(ej) .*...
        Np(ei - ej, Tph) .* (ones(length(e), 1) * de') / Tc^3;
    Gs(~isfinite(Gs)) = 0;

    Gs_in = Gs';
    Gs_out = sum(Gs, 2);
end

function Gr = Grecombination(e, Tph, Tc)
% Recombination matrix defined by Eq. (7) in J. M. Martinis et al.,
% Phys. Rev. Lett. 103, 097002 (2009), as well as Eqs. (C5) and (C6) in 
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
% These equations come from Eq. (8) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).
    [ej, ei] = meshgrid(e);

    Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) .* Np(ei + ej, Tph) / Tc^3;
end

function Gt = Gtrapping(e, Tph, Tc)
    % Simple model for the trapping matix.
    de = .5e-4; % in units of \Delta
    e_gap = de/2:de:1-de/2;
    [eg, ei] = meshgrid(e_gap, e);

    Gt = (ei - eg).^2 .* Np(ei - eg, Tph) / Tc^3;
    Gt = c * trapz(e_gap, Gt, 2);
end

function Gb = Gbreaking(e, de, Tph, Tc)
% Pair-breaking matrix, Eq. (27) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).

% Integration should be from \epsilon + \Delta to \Omega_D but
% because T_{ph} < T_c << \Omega_D we can truncate the integration earlier,
% e.g. at about 20 \Delta.
    dOmega = 2e-5; % in units of \Delta
    Omega = 2:dOmega:2+20*Tph;
    [Omega, ei] = meshgrid(Omega, e);

    Gb = (Omega.^2 ./ (exp(Omega / Tph) - 1)) .* ...
         (ei .* (Omega - ei) + 1) ./...
         (sqrt(ei.^2 - 1) .* sqrt((Omega - ei).^2 - 1)) / Tc^3;
    Gb(Omega < ei + 1) = 0;
    
    Z1_0 = 1.43; % See Table I in S. B. Kaplan et al.,
    % Phys. Rev. B 14, 4854 (1976). 
    Gb = (2 / Z1_0) * trapz(Omega, Gb, 2) .* de;
end

function R = Injection(e, rho_de, V, r)
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
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gb, Gt, R)
    if t > 0
        R = 0;
    end
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) +...
        Gb - Gt .* n + R;
end