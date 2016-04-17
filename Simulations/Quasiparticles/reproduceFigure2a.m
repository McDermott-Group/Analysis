function reproduceFigure2a
%reproduceFigure1 Reproduce Figure 2(a) from J. M. Martinis et al.,
%Phys. Rev. Lett. 103, 097002 (2009).

Tph = 0:.005:.180; % K
V = [2.8, 3.0]; % in units of \Delta
r = 2 * 1.7 * 1e-10 / (sqrt(8) - sqrt(2.8^2 - 1)); % in units of 1 / \tau_0
                                      %(assuming n_{qp} in units of n_{cp})
tspan = [-30000, 0]; % in units of tau0

n_qp_eq = nan(size(Tph));
for k = 1:length(Tph)
    [~, ~, ~, ~, ~, n_qp_T, r_qp] = noTrapping0DModel(r, V, Tph(k), tspan);
    n_qp_eq(k) = n_qp_T(end);
end

figure
plot(Tph, n_qp_eq / n_qp_eq(1), 'LineWidth', 3)
xlabel('Temperature T_p [K]', 'FontSize', 14)
ylabel('\Gamma_1(T) / \Gamma_1(0)', 'FontSize', 14)
title(['r_{qp} = ', num2str(r_qp, '%.2e'), ' / \tau_0'])
axis([0 .2 .95 1.15])
grid on

end