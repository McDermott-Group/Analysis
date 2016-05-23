function [t, e, n, f, n_qp, r_qp, P] = ...
    twoRegionTimeDomainModel(Tph, tspan, V, rqp, rph, c, vol, N)
%twoRegionTimeDomainModel Two-region, one that describes a normal metal-
% isolator-superconductor junction (NIS) and the other one - a resonator,
% quasi-0D model for computing the quasiparticle density evolution.  
% 
% [t, e, n, f, n_qp, r_qp, P] = 
%   twoRegionTimeDomainModel(Tph, tspan, V, rqp, rph, c, vol, N)
%   computes the quasiparticle dynamics using a two-region quasi-0D model.
%   The quasiparticles are directly injected into the NIS junction.
%   At every point in time the phonon density is computed. This phonon
%   density is used to calculate the quasiparticle injection rate assuming
%   the pair-breaking mechanism.
%
%   The input parameters:
%      Tph is phonon temperature in K,
%      tspan should be in form of [ti, tf] in units of \tau_0 where ti is
%      the inital time and tf is the final time of the integration range
%     (it is assumed quasiparticles are injected at t < 0 and relaxation
%      happens at t > 0),
%      V is the applied voltage in \Delta,
%      rqp is the quasiparticle injection rate in n_{cp}/\tau_0,
%      rph is the fraction of the phonons that go to the resonator,
%      dimensionless,
%      c is the trapping (capture) rate in n_{cp}/\tau_0,
%      vol is the injection volume in um^3,
%      N is the number of the energy bins (50-1000 or so).
%
%   The output parameters:
%      t is time in \tau_0,
%      e are energies of the bins in \Delta, a vector,
%      n are the quasiparticle densities in n_{cp}/\tau_0,
%      size(n) == [length(t), length(e)],
%      f are the occupational numbers, size(f) == size(n),
%      n_qp is the non-equilibrium quasiparticle density in quasiparticles
%      per um^3,
%      length(n_qp) == length(t),
%      r_qp is the total injection rate in n_{cp}/\tau_0, a single number,
%      P is the total injected power in W, a single number.

% Constants.
kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

delta = 0.18e-3; % eV (aluminum superconducting gap)

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta

Tc = 1 / 1.764; % \delta/(k_B * T_c) = 1.764 at T = 0 K, BCS result

tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 um^-3
           % C. Wang et al. Nature Comm. 5, 5836 (2014)

% Maximum energy.
max_e = max(4, max(V));

% Assign the quasiparicle energies to the bins. Non-uniform energy
% spacing is implemented. To get a spacing that is close to a uniform
% set alpha to a very small positive value.
alpha = 6;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

% Short-hand notation.
rho_de = rho(e) .* de;

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

Rdirect = DirectInjection(e, rho_de, V, rqp);

% Initial condition for the quasiparticle distribution.
% The first half of vector n0 describes the NIS junction and the second -
% the resonator.
n0 = zeros(size(e));
n0 = [n0; n0];

% Solve the ODE. 
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr,...
    Rdirect, e, de, rho_de, V, rph, Tc, Tph, c), tspan, n0, options);

n = n(:, size(n, 2)/2+1:end);
% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipartical density.
n_qp = 2 * ncp * sum(n, 2); % um^-3

% Bins indices the quasiparticles are injected to.
if length(V) == 2
    indices = V(1) < e & e < V(2);
else
    indices = e < V;
end

% Total injection rate.
r_qp = rqp * sum(rho_de(indices)); % n_{cp}/\tau_0

% Injection power.
P = rqp * sum(e(indices) .* rho_de(indices));
coeff = delta * ncp * vol / tau_0;
P = coeff * P; % W
end

function normalized_density = rho(e)
    normalized_density = e ./ sqrt(e.^2 - 1);
    normalized_density(e <= 1) = 0;
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

    Gs = (ei - ej).^2 .*...
        (1 - 1 ./ (ei .* ej)) .*...
        rho(ej) .* Np(ei - ej, Tph) .* (ones(length(e), 1) * de') / Tc^3;
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

    Gr = (ei + ej).^2 .*...
        (1 + 1 ./ (ei .* ej)) .*...
        Np(ei + ej, Tph) / Tc^3;
end

function Gtr = Gtrapping(e, Tph, Tc, c)
    % Simple model for the trapping matix.
    de = .5e-4; % in units of \Delta
    e_gap = de/2:de:1-de/2;
    [eg, ei] = meshgrid(e_gap, e);
    Gtr = (ei - eg).^2 .* Np(ei - eg, Tph) / Tc^3;
    Gtr = c * trapz(e_gap, Gtr, 2);
end

function R = DirectInjection(e, rho_de, V, r)
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation - at t > 0.
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

function [R, Omega1D, N_Omega] = ScatteringInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    persistent N Omega e_init e dN_Omega_no_f_inj dR_Omega_no_N_Omega
    % If the bias is too small there is no point in computing
    % the contribution due to the phonon scattering.
    if max(V) <= 3
        R = zeros(size(e_inj));
        Omega1D = [0, 10];
        N_Omega = [0, 0];
        return
    end
    % A few matrices are needed to be computed only once. Their value
    % will stay in the memory between the function calls since they are
    % intialized as "persistent".
    if isempty(N)
        N = 5 * length(e_inj);
        Omega = linspace(2, max(V) - 1, N);
        e_init = linspace(min(e_inj), max(e_inj), N)';
        [Omega, e] = meshgrid(Omega, e_init);
        % The following expression, multiplied by the quasiparticle
        % occupation numbers f_inj and integrated over the quasiparticle
        % "initial" energies e_init, gives phonon density per unit volume
        % at energy Omega.
        % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
        dN_Omega_no_f_inj = r * Omega.^2 .* rho(e - Omega) .* rho(e) .*...
                (1 - 1 ./ (e .* (e - Omega))) .* Np(Omega, Tph) / Tc^3;
        % The following expression, multiplied by the phonon density per 
        % unit volume and integrated over the phonon energies, gives
        % quasiparticle injection rates per unit energy per unit volume.
        % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
        dR_Omega_no_N_Omega = Omega.^2 .* (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1)) / Tc^3;
        dR_Omega_no_N_Omega(dR_Omega_no_N_Omega < 0 |...
                  ~isfinite(dR_Omega_no_N_Omega)) = 0;
    end
    dN_Omega = dN_Omega_no_f_inj .* f_inj(e);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = trapz(e_init, dN_Omega);

    dR = dR_Omega_no_N_Omega .* (ones(size(e_init)) * N_Omega);
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2);
    % Compute the quasipartcle injection per energy bin.
    R = interp1(e_init, R, e_inj) .* de_inj;
    Omega1D = Omega(1, :);
end

function [R, Omega1D, N_Omega] = RecombinationInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    persistent N Omega e_final e dN_Omega_no_f_inj dR_Omega_no_N_Omega
    % A few matrices are needed to be computed only once. Their value
    % will stay in the memory between the function calls since they are
    % intialized as "persistent".
    if isempty(N)
        N = 3 * length(e_inj);
        Omega = linspace(2, 2 * max(V), N);
        e_final = linspace(1, max(e_inj), N)';
        [Omega, e] = meshgrid(Omega, e_final);
        % The following expression, multiplied by the quasiparticle
        % occupation numbers f_inj at the recombininng energies and
        % integrated over the quasiparticle "final" energies e_final,
        % gives phonon density per unit volume at energy Omega.
        % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
        dN_Omega_no_f_inj = r * Omega.^2 .* rho(Omega - e) .* rho(e) .*...
                (1 + 1 ./ (e .* (Omega - e))) .* Np(Omega, Tph) / Tc^3;
        % The following expression, multiplied by the phonon density per 
        % unit volume and integrated over the phonon energies, gives
        % quasiparticle injection rates per unit energy per unit volume.
        % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
        dR_Omega_no_N_Omega = Omega.^2 .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1)) / Tc^3;
        dR_Omega_no_N_Omega(dR_Omega_no_N_Omega < 0 |...
                  ~isfinite(dR_Omega_no_N_Omega)) = 0;
    end
    dN_Omega = dN_Omega_no_f_inj .* f_inj(Omega - e) .* f_inj(e);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = trapz(e_final, dN_Omega);

    dR = dR_Omega_no_N_Omega .* (ones(size(e_final)) * N_Omega);
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2);
    % Compute the quasipartcle injection per energy bin.
    R = interp1(e_final, R, e_inj) .* de_inj;
    Omega1D = Omega(1, :);
end

function [R, Omega1D, N_Omega] = TrapInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    persistent N Omega e_init e dN_Omega_no_f_inj dR_Omega_no_N_Omega
    % If the bias is too small there is no point in computing
    % the contribution due to the phonon generation via trapping.
    if max(V) <= 2
        R = zeros(size(e_inj));
        Omega1D = [0, 10];
        N_Omega = [0, 0];
        return
    end
    % A few matrices are needed to be computed only once. Their value
    % will stay in the memory between the function calls since they are
    % intialized as "persistent".
    if isempty(N)
        N = 3 * length(e_inj);
        Omega = linspace(2, max(V), N);
        e_init = linspace(min(e_inj), max(e_inj), N)';
        [Omega, e] = meshgrid(Omega, e_init);
        % The following expression, multiplied by the quasiparticle
        % occupation numbers f_inj and integrated over the quasiparticle
        % "initial" energies e_init, gives phonon density per unit volume
        % at energy Omega.
        % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976),
        % assuming \Delta = 0.
        dN_Omega_no_f_inj = r * Omega.^2 .* rho(e) .* Np(Omega, Tph) / Tc^3;
        dN_Omega_no_f_inj(e - Omega >= 1 | e - Omega < 0) = 0;
        % The following expression, multiplied by the phonon density per 
        % unit volume and integrated over the phonon energies, gives
        % quasiparticle injection rates per unit energy per unit volume.
        % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
        dR_Omega_no_N_Omega = Omega.^2 .* (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1)) / Tc^3;
        dR_Omega_no_N_Omega(dR_Omega_no_N_Omega < 0 |...
                  ~isfinite(dR_Omega_no_N_Omega)) = 0;
    end
    dN_Omega = dN_Omega_no_f_inj .* f_inj(e);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = trapz(e_init, dN_Omega);

    dR = dR_Omega_no_N_Omega .* (ones(size(e_init)) * N_Omega);
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2);
    % Compute the quasipartcle injection per energy bin.
    R = interp1(e_init, R, e_inj) .* de_inj;
    Omega1D = Omega(1, :);
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr, R_direct,...
        e, de, rho_de, V, rph, Tc, Tph, c)
    % n is in units n_{cp}, rates are in units n_{cp}/\tau_0.
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation - at t > 0.
    if t > 0
        R_direct = 0;
    end
    non_positive_n = n <= 0;
    n(non_positive_n) = 0;

    % The first half of vector n0 describes the NIS junction and
    % the second - the resonator.
    half = length(n) / 2;
    n_cnt = n(1:half);
    n_res = n(half+1:end);

    f_inj = n_cnt ./ rho_de;
    f_inj(f_inj < 0) = 0;
    f_inj = @(e_inj) interp1(e, f_inj, e_inj);

    Rph_rec = RecombinationInjection(e, de, f_inj, V, rph, Tc, Tph);
    Rph_sct = ScatteringInjection(e, de, f_inj, V, rph, Tc, Tph);
    Rph_trp = TrapInjection(e, de, f_inj, V, c * rph, Tc, Tph);

    ndot_cnt = Gs_in * n_cnt - Gs_out .* n_cnt - 2 * n_cnt .* (Gr * n_cnt) -...
        Gtr .* n_cnt + R_direct;
    ndot_res = Gs_in * n_res - Gs_out .* n_res - 2 * n_res .* (Gr * n_res) -...
        Gtr .* n_res + Rph_rec + Rph_trp + Rph_sct;
    ndot = [ndot_cnt; ndot_res];
    ndot(ndot < 0 & non_positive_n) = 0;
end