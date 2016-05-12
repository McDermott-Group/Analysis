function [t, e, n, f, n_qp, r_qp, P] = ...
    twoStageQuasi0DModel(Tph, tspan, V, rqp, rph, c, plot_flag)
%twoStageQuasi0DModel(Simple trapping, quasi-0D model.
% [t, e, n, f, n_qp, r_qp, P] = 
%   phononMediatedQuasi0DModel(Tph, tspan, V, r, c, plot_flag)
%   computes the quasiparticle dynamics in a simple model without any
%   quasiparticle traps.
%   The input parameters:
%      Tph is phonon temperature in K,
%      tspan should be in form of [ti, tf] in units of \tau_0 where ti is
%      the inital time and tf is the final time of the integration (it is
%      assumed quasiparticles are injected at t < 0 and relaxation happens
%      at t > 0),
%      V is the applied voltage in \Delta,
%      r_inj is the quasiparticle injection rate in 1/\tau_0,
%      assuming n and n_qp in units of n_{cp},
%      r_ph is the injection rate in 1/\tau_0, assuming n and n_qp in units
%      of n_{cp},
%      c is the trapping (capture) rate in 1/\tau_0, assuming n and n_qp
%      in units of n_{cp},
%      plot_flag is an optional parameter that controls plot creation.
%   The output parameters:
%      t is time in \tau_0,
%      e are energies of the bins in \Delta,
%      n are the quasiparticle densities in n_{cp}, size(n) == [length(t),
%      length(e)],
%      f are the occupational numbers, size(f) == size(n),
%      n_qp is the non-equilibrium quasiparticle density,
%      length(n_qp) == length(t),
%      r_qp is the total injection rate in 1 / \tau_0, assuming n and n_qp
%      in units n_{cp}, single number,
%      P is the total injected power in W, single number.

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

Volume = 1e5; % um^3
tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)
                     
% Debye temperature.
% TD = 428; % K
% TD = kB * TD/ delta; % in units of \Delta

% Number of energy bins.
N = 500;
% Maximum energy.
max_e = 2 * max(V) - 1;

% Assign the quasiparicle energies to the bins. Non-uniform energy
% separation is implemented. For a close to uniform distribution set alpha
% to a small positive value.
%e = linspace(1, max_e, N + 1)';

alpha = 5;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

% Short-hand notation.
rho_de = rho(e) .* de;

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

Rqp = DirectInjection(e, rho_de, V, rqp);

% Initial condition n0.
% f0 = 1 ./ (exp(e / Tph) + 1);
% n0_T = rho_de .* f0;
n0 = zeros(size(e));

% Solve the ODE. 
options = odeset('AbsTol', 1e-20);
[~, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rqp), tspan, n0, options);

n_inj = n(end, :);

% Occupational numbers.
f_inj = n_inj ./ rho_de';
f_inj(f_inj < 0) = 0;
f_inj = @(e_inj) interp1(e, f_inj, e_inj, 'spline', 'extrap');

Rph_rec = RecombinationInjection(e, de, f_inj, V, rph, Tc, Tph);
Rph_sct = ScatteringInjection(e, de, f_inj, V, rph, Tc, Tph);

Gtr = Gtrapping(e, Tph, Tc, c);

% Generate a plot.
if plot_flag
    f0 = 1 ./ (exp(e / Tph) + 1);
    n0 = rho_de .* f0;
    figure
    plot(e, Gs_in * n0 ./ de, e,  Gs_out .* n0 ./ de,...
         e, n0 .* ( Gr * n0) ./ de,...
         e, Gtr .* n0 ./ de, e, Rph_sct ./ de, e, Rph_rec ./ de, 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Rate per Unit Energy (1/\tau_0\Delta)', 'FontSize', 14)
    legend({'G^{in}_{scattering}', 'G^{out}_{scattering}',...
            'G_{recombination}', 'G_{trapping}', 'R_{ph. sct.}', 'R_{ph. rec.}'})
    title('Absolute Term Strengths at Thermal Equilibrium')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    
    figure
    plot(e, Rph_sct, e, Rph_rec, 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Injection Rate per Unit Energy (1/\tau_0\Delta)', 'FontSize', 14)
    legend({'R_{ph. sct.}', 'R_{ph. rec.}'})
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
end

options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rph_sct + Rph_rec), tspan, n0, options);

% figure
% loglog(e, n_inj, e, n(end, :))
% legend('Direct Injection', 'Solution')

% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipatical density.
n_qp = 2 * sum(n, 2);

% Bins indices the quasiparticles are injected to.
if length(V) == 2
    indices = V(1) < e & e < V(2);
else
    indices = e < V;
end

% Total injection rate.
r_qp = rqp * sum(rho_de(indices));

% Injection power.
P = rqp * sum(e(indices) .* rho_de(indices));
coeff = delta * ncp * Volume / tau_0;
P = coeff * P; % in W
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

function [R, Omega, N_Omega] = ScatteringInjection(e_inj, de_inj, f_inj, V, r, Tc, Tph)
    N = 5 * length(e_inj);
    Omega = linspace(0, max(V) - 1, N);
    e_final = linspace(min(e_inj), max(e_inj), N)';
    [Omega, e] = meshgrid(Omega, e_final);
    dN_Omega = Omega.^2 .* rho(e - Omega) .* rho(e) .* f_inj(e) .*...
        (1 - 1 ./ (e .* (e - Omega))) .* Np(Omega, Tph);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = r * trapz(e_final, dN_Omega) / Tc^3;
    
%     figure
%     plot(Omega(1, :), N_Omega, 'LineWidth', 2)
%     xlabel('Phonon Energy \Omega (\Delta)', 'FontSize', 14)
%     ylabel('N_{gen}(\Omega) (n_{cp} / \Delta)', 'FontSize', 14)
%     title('Phonon Spectrum (Scattering)')
%     set(gca, 'yscale', 'Log')
%     axis tight
%     grid on

    dR = Omega.^2 .* (ones(size(e_final)) * N_Omega) .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1));
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2) / Tc^3;
    R = interp1(e_final, R, e_inj) .* de_inj;
    Omega = Omega(1, :);
end

function [R, Omega, N_Omega] = RecombinationInjection(e_inj, de_inj, f_inj, V, r, Tc, Tph)
    N = 5 * length(e_inj);
    Omega = linspace(2, 2 * max(V), N);
    e_final = linspace(min(e_inj), max(e_inj), N)';
    [Omega, e] = meshgrid(Omega, e_final);
    dN_Omega = Omega.^2 .* rho(Omega - e) .* rho(e) .*...
        f_inj(Omega - e) .* f_inj(e) .*...
        (1 + 1 ./ (e .* (Omega - e))) .* Np(Omega, Tph);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = r * trapz(e_final, dN_Omega) / Tc^3;

%     figure
%     plot(Omega(1, :), N_Omega, 'LineWidth', 2)
%     xlabel('Phonon Energy \Omega (\Delta)', 'FontSize', 14)
%     ylabel('N_{gen}(\Omega) (n_{cp} / \Delta)', 'FontSize', 14)
%     title('Phonon Spectrum (Recombination)')
%     set(gca, 'yscale', 'Log')
%     axis tight
%     grid on

    dR = Omega.^2 .* (ones(size(e_final)) * N_Omega) .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1));
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2) / Tc^3;
    R = interp1(e_final, R, e_inj) .* de_inj;
    Omega = Omega(1, :);
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr, R)
    if t > 0
        R = 0;
    end
    n(n < 0) = 0;
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) - Gtr .* n + R;
end