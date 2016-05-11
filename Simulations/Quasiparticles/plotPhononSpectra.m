function plotPhononSpectra
%plotPhononSpectra Plot phonon spectra.

r_direct = .09e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 1; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 0; % trapping rate in units of 1 / \tau_0

Tph = 0.051; % K
tspan = [-500, 0]; % in units of \tau_0

V = 1.1:1:10;

Ptot = NaN(size(V));
Psct = NaN(size(V));
Prec = NaN(size(V));
Ptotsct = NaN(size(V));
for k = 1:length(V)
    if V(k) > 1
        [~, ~, ~, ~, ~, ~, Ptot(k), Psct(k), Prec(k), Ptotsct(k)] = phononSpectra(Tph, tspan,...
            V(k), r_direct, r_phonon, c);
    end
    k
end

figure
plot(V, Ptotsct ./ Ptot, V, Psct ./ Ptot, V, Prec ./ Ptot, V, (Psct + Prec) ./ Ptot,...
    V, (Ptotsct + Prec) ./ Ptot, 'LineWidth', 2)
xlabel('Voltage (\Delta)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'scattering above 2\Delta', 'recombination',...
    'scattering above 2\Delta + recombination', 'scattering + recombination')
axis tight
grid on

figure
semilogx(Ptot, Ptotsct ./ Ptot, Ptot, Psct ./ Ptot, Ptot, Prec ./ Ptot, Ptot,...
    (Psct + Prec) ./ Ptot, Ptot, (Ptotsct + Prec) ./ Ptot,'LineWidth', 2)
xlabel('P_{total} (W)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'scattering above 2\Delta', 'recombination',...
    'scattering above 2\Delta + recombination', 'scattering + recombination')
axis tight
grid on

end