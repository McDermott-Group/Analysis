function [t, e, n, f, n_qp, r_qp, Ptot, Psct, Prec, Psct2D,...
    P_sct, P_rec] = ...
    phononSpectra(Tph, tspan, V, rqp, rph, c, vol, N)
%phononSpectra(Simple trapping, quasi-0D model.
% [t, e, n, f, n_qp, r_qp, P] = 
%   phononMediatedQuasi0DModel(Tph, tspan, V, r, c)
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
%      vol is the volume in um^3,
%      N is the umber of energy bins (50-1000 or so),
%      plot_flag is an optional parameter that controls plot creation.
%   The output parameters:
%      t is time in \tau_0,
%      e are energies of the bins in \Delta,
%      n are the quasiparticle densities in  in um^-3,, size(n) == [length(t),
%      length(e)],
%      f are the occupational numbers, size(f) == size(n),
%      n_qp is the non-equilibrium quasiparticle density,
%      length(n_qp) == length(t),
%      r_qp is the total injection rate in 1 / \tau_0, single number,
%      P is the total injected power in \Delta / \tau_0, single number.

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

tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)

power_calib = delta * ncp * vol / tau_0;
                     
% Debye temperature.
% TD = 428; % K
% TD = kB * TD/ delta; % in units of \Delta

% Maximum energy.
max_e = 2 * max(V) - 1;

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

% Initial condition n0.
% f0 = 1 ./ (exp(e / Tph) + 1);
% n0_T = rho_de .* f0;
n0 = zeros(size(e));

% Solve the ODE. 
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rqp), tspan, n0, options);

n_inj = n(end, :);

% Occupational numbers.
f_inj = n_inj ./ rho_de';
f_inj(f_inj < 0) = 0;
f_inj = @(e_inj) interp1(e, f_inj, e_inj);

[~, Omega_rec, N_Omega_rec] = RecombinationInjection(e, de, f_inj, V, rph, Tc, Tph);
[~, Omega_sct, N_Omega_sct] = ScatteringInjection(e, de, f_inj, V, rph, Tc, Tph);

% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipatical density.
n_qp = 2 * ncp * sum(n, 2);

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
Ptot = power_calib * P; % in W
Psct = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);
N_Omega_sct = N_Omega_sct(Omega_sct > 2);
Omega_sct = Omega_sct(Omega_sct > 2);
Psct2D = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);
Prec = power_calib * trapz(Omega_rec, N_Omega_rec .* Omega_rec);

P_sct_in = power_calib * sum(e .* (Gs_in * n_inj'));
P_sct_out = power_calib * sum(e .* Gs_out .* n_inj');
P_sct = P_sct_out - P_sct_in;
P_rec = power_calib * sum(2 * e .* n_inj' .* (Gr * n_inj'));

Pt_sct = NaN(size(n, 1), 1);
Pt_sct2D = NaN(size(n, 1), 1);
Pt_rec = NaN(size(n, 1), 1);
f(f < 0) = 0;
for k = 1:size(n, 1)
    f_inj = f(k, :);
    f_inj = @(e_inj) interp1(e, f_inj, e_inj);
    [~, Omega_rec, N_Omega_rec] = RecombinationInjection(e, de, f_inj, V, rph, Tc, Tph);
    [~, Omega_sct, N_Omega_sct] = ScatteringInjection(e, de, f_inj, V, rph, Tc, Tph);
    Pt_sct(k) = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);
    N_Omega_sct = N_Omega_sct(Omega_sct > 2);
    Omega_sct = Omega_sct(Omega_sct > 2);
    Pt_sct2D(k) = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);
    Pt_rec(k) = power_calib * trapz(Omega_rec, N_Omega_rec .* Omega_rec);
end

figure
t = t - min(t);
plot(t, Pt_sct ./ Ptot, t, Pt_rec ./ Ptot, t, (Pt_sct + Pt_rec) ./ Ptot,...
     t, Pt_sct2D ./ Ptot, t, (Pt_sct2D + Pt_rec) ./ Ptot, 'LineWidth', 2)
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination',...
       'scattering + recombination', 'scattering above 2\Delta',...
       'scattering above 2\Delta + recombination')
title(['V = ',  num2str(V), ' \Delta/e']) 
axis tight
grid on

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
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).c
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
    non_positive_n = n <= 0;
    n(non_positive_n) = 0;
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) - Gtr .* n + R;
    ndot(ndot < 0 & non_positive_n) = 0;
end