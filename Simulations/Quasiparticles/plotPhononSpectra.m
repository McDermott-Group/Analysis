function plotPhononSpectra
%plotPhononSpectra Plot phonon spectra.

r_direct = .09e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon = 1; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c = 0; % trapping rate in units of 1 / \tau_0

Tph = 0.051; % K
tspan = [-500, 0]; % in units of \tau_0

V = 1.1:2:10;

Ptot = NaN(size(V));
Psct2D = NaN(size(V));
Prec = NaN(size(V));
Psct = NaN(size(V));
P_sct = NaN(size(V));
P_rec = NaN(size(V));
for k = 1:length(V)
    if V(k) > 1
        [~, ~, ~, ~, ~, ~, Ptot(k), Psct(k), Prec(k), Psct2D(k),...
            P_sct(k), P_rec(k)] = phononSpectra(Tph, tspan,...
            V(k), r_direct, r_phonon, c);
    end
    k
end

figure
plot(V, Psct ./ Ptot, V, Prec ./ Ptot, V, (Psct + Prec) ./ Ptot,...
     V, Psct2D ./ Ptot, V, (Psct2D + Prec) ./ Ptot, 'LineWidth', 2)
xlabel('Voltage (\Delta)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination',...
       'scattering + recombination', 'scattering above 2\Delta',...
       'scattering above 2\Delta + recombination')
axis tight
grid on

figure
plot(V, P_sct ./ Ptot, V, P_rec ./ Ptot,...
    V, (P_sct + P_rec) ./ Ptot, 'LineWidth', 2)
xlabel('Voltage (\Delta)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination',...
    'scattering + recombination')
axis tight
grid on

figure
semilogx(Ptot, Psct ./ Ptot, Ptot, Prec ./ Ptot, Ptot, (Psct + Prec) ./ Ptot,...
     Ptot, Psct2D ./ Ptot, Ptot, (Psct2D + Prec) ./ Ptot, 'LineWidth', 2)
xlabel('P_{total} (W)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination',...
       'scattering + recombination', 'scattering above 2\Delta',...
       'scattering above 2\Delta + recombination')
axis tight
grid on

figure
semilogx(Ptot, P_sct ./ Ptot, Ptot, P_rec ./ Ptot, Ptot,...
    (P_sct + P_rec) ./ Ptot, 'LineWidth', 2)
xlabel('P_{total} (W)', 'FontSize', 14)
ylabel('Normalized Power P / P_{total}', 'FontSize', 14)
legend('scattering', 'recombination',...
    'scattering + recombination')
axis tight
grid on

end