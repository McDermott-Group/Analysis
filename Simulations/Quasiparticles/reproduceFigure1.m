function reproduceFigure1
%reproduceFigure1 Reproduce Figure 1 from J. M. Martinis et al.,
%Phys. Rev. Lett. 103, 097002 (2009).

r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1)); % in units of 1 / \tau_0
                                      %(assuming n_{qp} in units of n_{cp})
V = [2.8, 3.0]; % in units of \Delta
Tph = [0.140, 0.070, 0]; % K
tspan = [-30000, 0]; % in units of \tau_0

figure
hold on
n_qp_eq = nan(size(Tph));
ax = gca;
ax.ColorOrderIndex = 2; 
for k = 1:length(Tph)
    [~, e, ~, f, n_qp, ~, r_qp] = noTrapping0DModel(r, V, Tph(k), tspan);
    n_qp_eq(k) = n_qp(end);
    plot(e, f(end, :), 'LineWidth', 3)
end

% Thermal distribution.
kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

delta = 0.18e-3; % eV (aluminum superconducting gap)

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = 0.140; % K
Tph = kB * Tph / delta; % in units of \Delta

for k = 1:length(Tph)
    fT = 1 ./ (exp(e / Tph(k)) - 1);
    plot(e, fT, 'LineWidth', 2)
end

xlabel('Quasiparticle energy E/\Delta', 'FontSize', 14)
ylabel('State occupation f(E)', 'FontSize', 14)
title({['n_{qp} / 2D(E_F)\Delta = ', num2str(n_qp_eq(1)/2, '%.2e'),...
        '; r_{qp} = ', num2str(r_qp, '%.2e'), ' / \tau_0']})
legend({'T_p = 140 mK', '          70 mK', '            0',...
        'f_T [140 mK]'})
set(gca, 'yscale', 'Log')
axis([1 2.4 1e-12 1.01e-5])
grid on

end