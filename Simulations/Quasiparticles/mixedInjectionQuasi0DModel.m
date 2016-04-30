function [t, e, n, f, n_qp, r_qp, P] = ...
    mixedInjectionQuasi0DModel(Tph, tspan, V, rqp, rph, c, plot_flag)
%mixedInjectionQuasi0DModel Simple trapping, quasi-0D model.
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
%      P is the total injected power in \Delta / \tau_0, assuming n and
%      n_qp in units n_{cp}, single number.

if ~exist('plot_flag', 'var')
    plot_flag = false;
end

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
N = 100;
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

Gr = Grecombination(e, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

Rqp = DirectInjection(e, rho_de, V, rqp);
Rph = PhononMediatedInjection(e, de, V, rph, Tc);

% Generate a plot.
if plot_flag
    f0 = 1 ./ (exp(e / Tph) + 1);
    n0 = rho_de .* f0;
    figure
    plot(e, Gs_in * n0 ./ de, e,  Gs_out .* n0 ./ de,...
         e, n0 .* ( Gr * n0) ./ de,...
         e, Gtr .* n0 ./ de, e, Rqp ./ de, e, Rph ./ de, 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Rate per Unit Energy (1/\tau_0\Delta)', 'FontSize', 14)
    legend({'G^{in}_{scattering}', 'G^{out}_{scattering}',...
            'G_{recombination}', 'G_{trapping}', 'R_{direct}', 'R_phonon'})
    title('Absolute Term Strengths at Thermal Equilibrium')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    
    figure
    plot(e, Rqp, e, Rph, 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Injection Rate per Unit Energy (1/\tau_0\Delta)', 'FontSize', 14)
    legend({'R_{direct}', 'R_{phonon}'})
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    return
end

% Initial condition n0.
% f0 = 1 ./ (exp(e / Tph) + 1);
% n0_T = rho_de .* f0;
n0 = zeros(size(e));

% Solve the ODE.
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rqp + Rph), tspan, n0, options);

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

function R = PhononMediatedInjection(e_inj, de_inj, V, r, Tc)
    N = 5e3;
    epsilon = 1e-4;
    Omega = linspace(0, V, N);
    N_Omega = zeros(size(Omega));
    for k = 1:length(Omega)
%         if Omega(k) < V - 1
%             e = linspace(1 + epsilon + Omega(k), V, N);
%             dN_Omega = Omega(k)^2 * rho(e - Omega(k)) .* rho(e) .*...
%                 (1 - 1 ./ (e .* (Omega(k) - e)));
%             N_Omega(k) = r * trapz(e, dN_Omega) / Tc^3;
%         end

        le = max(1 + epsilon, Omega(k));
        he = min(V, Omega(k) + 1);
        if le < he
            e = linspace(le, he, N);
            dN_Omega = Omega(k)^2 * rho(e);
            N_Omega(k) = r * trapz(e, dN_Omega) / Tc^3;
        end
    end
    R = zeros(size(e_inj));
    for k = 1:length(e_inj)
        if e_inj(k) <= V - 1
           indices = Omega > (e_inj(k) + 1 + epsilon);
           if length(Omega(indices)) > 1
                dR = Omega(indices).^2 .* N_Omega(indices) .*...
                    (e_inj(k) * (Omega(indices) - e_inj(k)) + 1) ./...
                    (sqrt(e_inj(k)^2 - 1) .*...
                    sqrt((Omega(indices) - e_inj(k)).^2 - 1));
                R(k) = trapz(Omega(indices), dR);
            end
        end
    end
    R = R .* de_inj / Tc^3;
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr, R)
    if t > 0
        R = 0;
    end
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) - Gtr .* n + R;
end