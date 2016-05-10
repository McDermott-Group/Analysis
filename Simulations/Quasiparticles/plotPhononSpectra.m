function plotPhononSpectra
%plotPhononSpectra Plot phonon spectra.

r_direct = .09e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 1; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 0; % trapping rate in units of 1 / \tau_0

delta = 0.18e-3; % eV (aluminum superconducting gap)
e = 1.60217662e-19; % C

ncp = 4e6; % n_{cp} for aluminum is 4e-6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)

Tph = 0.051; % K
tspan = [-500, 0]; % in units of \tau_0

V = 1:.25:10;

Ptot = NaN(size(V));
Psct = NaN(size(V));
Prec = NaN(size(V));
for k = 1:length(V)
    if V(k) > 1
        [~, ~, ~, ~, ~, ~, Ptot(k), Psct(k), Prec(k)] = phononSpectra(Tph, tspan,...
            V(k), r_direct, r_phonon, c);
    end
    k
end

figure
plot(V, Psct ./ Ptot, V, Prec ./ Ptot, 'LineWidth', 2)
xlabel('Voltage (\Delta)', 'FontSize', 14)
ylabel('Normailized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination')
% set(gca, 'yscale', 'Log')
axis tight
grid on

figure
semilogx(Ptot, Psct ./ Ptot, Ptot, Prec ./ Ptot, 'LineWidth', 2)
xlabel('P_{total} (W)', 'FontSize', 14)
ylabel('Normailized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination')
% set(gca, 'yscale', 'Log')
axis tight
grid on

end