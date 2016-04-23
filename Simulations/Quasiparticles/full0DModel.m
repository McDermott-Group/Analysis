function [t, e, n, f, n_qp, n_qp_T, r_qp, P] = ...
    full0DModel(Tph, tspan, V, r, c, plot_flag)
%full0DModel Full 0D model.
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
%      plot_flag is an optional parameter that controls plot creation.
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
N = 200;
% Maximum energy.
max_e = max(V) + 1;

% Assign the quasiparicle energies to the bins. Non-uniform energy
% separation is implemented. For a close to uniform distribution set alpha
% to a small positive value.
alpha = 5;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

% Short-hand notation.
rho_de = rho(e) .* de;

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

[Gbs_in, Gbs_out] = Gbranch_scattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gbr = Gbranch_recombination(e, Tph, Tc);

Gpb = Gpair_breaking(e, de, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

R = Injection(e, rho_de, V, r);

% Initial condition n0.
f0 = 1 ./ (exp(e / Tph) + 1);
n0 = rho_de .* f0;

% Generate a plot.
if exist('plot_flag', 'var') && plot_flag
    figure
    Total =  Gs_in * n0 - Gs_out .* n0 + Gbs_in * n0 - Gbs_out .* n0 -...
         n0 .* ( Gr * n0) - n0 .* (Gbr * n0) +...
         Gpb - Gtr .* n0 + R;
    plot(e,  Gs_in * n0, e,  Gs_out .* n0,...
         e, Gbs_in * n0, e, Gbs_out .* n0,...
         e, n0 .* ( Gr * n0),...
         e, n0 .* (Gbr * n0),...
         e, Gpb, e, Gtr .* n0, e, R, e, abs(Total), 'k', 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Rate (1/\tau_0)', 'FontSize', 14)
    legend({'G^{in}_{intra-scattering}', 'G^{out}_{intra-scattering}',...
            'G^{in}_{inter-scattering}', 'G^{out}_{inter-scattering}',...
            'G_{intra-recombination}', 'G_{inter-recombination}',...
            'G_{pair-breaking}', 'G_{trapping}', 'R_{injection}', 'Total'})
    title('Absolute Term Strengths at Thermal Equilibrium')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    
    figure
    plot(e, Total ./ (n0 .* ( Gr * n0)),...
         e, Total ./ (n0 .* (Gbr * n0)),...
         e, Total ./ Gpb, e, Total ./ (Gtr .* n0), 'k', 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Rate_{total} / Rate_{process}', 'FontSize', 14)
    legend({'G_{intra-recombination}', 'G_{inter-recombination}',...
            'G_{pair-breaking}', 'G_{trapping}',})
    title('Relative Term Strengths at Thermal Equilibrium')
    axis tight
    ylim([-3, 3])
    grid on
    return
end

% Generate branching matrices (and vectors).
n0 = [n0; n0];
rho_de = [rho_de; rho_de];

Gs_in = [ Gs_in, Gbs_in;...
         Gbs_in,  Gs_in];
     
Gs_out = [ Gs_out, Gbs_out;...
          Gbs_out,  Gs_out];

Gr = [ Gr, Gbr;
      Gbr,  Gr];
  
Gpb = [Gpb; Gpb];
Gtr = [Gtr; Gtr];
R = [R; zeros(size(R))];

% Solve the ODE.
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gpb, Gtr, R), tspan, n0, options);

% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipatical density.
n_qp = sum(n, 2);

% Thermal quasipatical density.
n_qp_T = sum(rho_de ./ (exp(e / Tph) + 1));

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

    Gs = (ei - ej).^2 .*...
        (1 - 1 ./ (ei .* ej)) .*...
        rho(ej) .* Np(ei - ej, Tph) .* (ones(length(e), 1) * de') / Tc^3;
    Gs(~isfinite(Gs)) = 0;

    Gs_in = Gs';
    Gs_out = sum(Gs, 2);
end

function [Gbs_in, Gbs_out] = Gbranch_scattering(e, de, Tph, Tc)
% Inter-branch scattering matrix, Eq. (21) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976), the first two lines.
    [ej, ei] = meshgrid(e);

    Gbs = .5 * (ei - ej).^2 .*...
        (1 - (1 + sqrt((ei.^2 - 1) .* (ej.^2 - 1))) ./ (ei .* ej)) .*...
        rho(ej) .* Np(ei - ej, Tph) .* (ones(length(e), 1) * de') / Tc^3;
    Gbs(~isfinite(Gbs)) = 0;

    Gbs_in = Gbs';
    Gbs_out = sum(Gbs, 2);
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

function Gbr = Gbranch_recombination(e, Tph, Tc)
% Inter-branch recombination matrix, Eq. (21) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976), the last line.
    [ej, ei] = meshgrid(e);

    Gbr = .5 * (ei + ej).^2 .*...
        (1 + (1 - sqrt(ei.^2 - 1) .* sqrt(ej.^2 - 1)) ./ (ei .* ej)) .*...
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

function Gpb = Gpair_breaking(e, de, Tph, Tc)
% Pair-breaking matrix, Eq. (27) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).

% Integration should be from \epsilon + \Delta to \Omega_D but
% because T_{ph} < T_c << \Omega_D we can truncate the integration earlier,
% e.g. at about 50 \Delta.
%     de_phonon = .5e-4; % in units of \Delta
%     e_phonon = 2:de_phonon:2+100*Tph;
%     [Omega, ei] = meshgrid(e_phonon, e);
%     Gpb = (Omega.^2 ./ (exp(Omega / Tph) - 1)) .* ...
%          (ei .* (Omega - ei) + 1) ./...
%          (sqrt(ei.^2 - 1) .* sqrt((Omega - ei).^2 - 1)) / Tc^3;
%     Gpb(Omega < ei + 1) = 0;
%     Z1_0 = 1.43; % See Table I in S. B. Kaplan et al.,
%     % Phys. Rev. B 14, 4854 (1976). 
%     Gpb = (2 / Z1_0) * trapz(e_phonon, Gpb, 2) .* de;
    Gpb = zeros(size(e));
    alpha = 4;
    N = 100000;
    for k = 1:length(e)
        % Omega = linspace(e(k) + 1.00001, 90 * Tph, 1e5);
        Omega = e(k) + 1 + (100 * Tph * sinh(alpha * (1:N) / N) / sinh(alpha))';
        g_pb = (Omega.^2 ./ (exp(Omega / Tph) - 1)) .* ...
         (e(k) * (Omega - e(k)) + 1) ./...
         (sqrt(e(k)^2 - 1) .* sqrt((Omega - e(k)).^2 - 1));
        Gpb(k) = trapz(Omega, g_pb);
    end
    Z1_0 = 1.43; % See Table I in S. B. Kaplan et al.,
    % Phys. Rev. B 14, 4854 (1976). 
    Gpb = (2 / Z1_0) * Gpb .* de / Tc^3;

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

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gpb, Gtr, R)
    if t > 0
        R = 0;
    end
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) +...
        Gpb - Gtr .* n + R;
end