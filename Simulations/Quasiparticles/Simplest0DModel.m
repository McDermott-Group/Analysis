function [t, e, n, f, n_qp, r_qp, P] = simplest0DModel(r, V, Tph, tspan)

kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

Tel = 0.050; % mK
% Tc = 1.2; % K (aluminum critical temperature)

delta = 0.18e-3; % eV (aluminum superconducting gap) 

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta
Tel = kB * Tel / delta;

% Tc = kB * Tc / delta;
Tc = 1 / 1.764; % \delta/(K_B * T_c) = 1.764 at T = 0 K BCS result

N = 50;
max_e = V + 3;
de = (max_e - 1) / N;
e = (linspace(1 + de / 2, max_e - de / 2, N))';
[ej, ei] = meshgrid(e);

gs = (ei - ej).^2 .* (1 - 1 ./ (ei .* ej)) .* rho(ej) ./...
    abs(exp(-(ei - ej) / Tph) - 1) / Tc^3;
gs(~isfinite(gs)) = 0;

gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) .* rho(ej) ./...
    abs(exp(-(ei + ej) / Tph) - 1) / Tc^3;

rho_e = rho(e);

% It is assmumed that the injection is happening at t < 0 and
% the relaxation at t > 0.
n0 = rho_e ./ (1 + exp(e / Tel));

% Injection voltage.
R = zeros(size(ej));
R(ej <= V) = r;

[t, n] = ode45(@(t, n) quasiparticleODE(t, n, rho_e, gs, gr, R), tspan, n0);

n_qp = 2 * sum(n, 2) * de;
f = 2 * n ./ (ones(length(t), 1) * rho_e');

r_qp = r * sum(rho_e(e < V)) * de;
P = r * sum(e(e < V) .* rho_e(e < V)) * de;
end

function normalized_density = rho(e)
    normalized_density = e ./ sqrt(e.^2 - 1);
end

function ndot = quasiparticleODE(t, n, rho, gs, gr, R)
    f = ones(length(n), 1) * (n ./ rho)';
    Gs = gs .* (1 - f);
    Gr = gr .* f;
    if t >= 0
        R = 0;
    end
    ndot =  Gs' * n - sum(Gs, 2) .* n - 2 * n .* (Gr * n) + R * rho;
end