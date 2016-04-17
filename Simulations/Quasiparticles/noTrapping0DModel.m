function [t, e, n, f, n_qp, n_qp_T, r_qp, P] = noTrapping0DModel(r, V, Tph, tspan)

kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

Tqp_T = Tph; % K
% Tc = 1.2; % K (aluminum critical temperature)

delta = 0.18e-3; % eV (aluminum superconducting gap) 

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta
Tqp_T = kB * Tqp_T / delta;

% Tc = kB * Tc / delta;
Tc = 1 / 1.764; % \delta/(K_B * T_c) = 1.764 at T = 0 K BCS result

N = 1000;
max_e = max(V) + 1;

alpha = 2;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

[ej, ei] = meshgrid(e);

Gs = (ei - ej).^2 .* (1 - 1 ./ (ei .* ej)) .* rho(ej) ./...
    abs(exp(-(ei - ej) / Tph) - 1) / Tc^3 .*...
    (ones(length(e), 1) * de');
Gs(~isfinite(Gs)) = 0;

Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) ./...
    abs(exp(-(ei + ej) / Tph) - 1) / Tc^3;

rho_de =  rho(e) .* de;

% It is assmumed that the injection is happening at t < 0 and
% the relaxation at t > 0.
n_T = rho_de ./ (1 + exp(e / Tqp_T));
% n0 = zeros(size(rho_de));
n0 = n_T;

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

options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n, rho_de, Gs, Gr, R, n_T),...
    tspan, n0, options);

n_qp = 2 * sum(n - ones(length(t), 1) * n_T', 2);
n_qp_T = 2 * sum(n, 2);
f = n ./ (ones(length(t), 1) * rho_de');

if length(V) == 2
    indices = V(1) < e & e < V(2);
else
    indices = e < V;
end

r_qp = r * sum(rho_de(indices));
P = r * sum(e(indices) .* rho_de(indices));

end

function normalized_density = rho(e)
    normalized_density = e ./ sqrt(e.^2 - 1);
end

function ndot = quasiparticleODE(t, n, rho_de, Gs, Gr, R, n_T)
    Gs = Gs .* (1 - ones(length(n), 1) * (n ./ rho_de)');
    if t > 0
        R = 0;
    end
    n = n - n_T;
    ndot = Gs' * n - sum(Gs, 2) .* n - 2 * n .* (Gr * n) + R;
end