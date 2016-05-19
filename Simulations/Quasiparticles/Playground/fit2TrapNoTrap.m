function fit2TrapNoTrap
%fit2TrapNoTrap Fitting to TrapNoTrap dataset.

r_direct_no_tr = 3e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_direct_no_tr = 0; % trapping rate in units of 1 / \tau_0

r_direct_tr = 3e-5; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_direct_tr = .08; % trapping rate in units of 1 / \tau_0

r_phonon_no_tr = 2.5e-9; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_phonon_no_tr = 0; % trapping rate in units of 1 / \tau_0

r_phonon_tr = 2.5e-9; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_phonon_tr = .2; % trapping rate in units of 1 / \tau_0

delta = 0.18e-3; % eV (aluminum superconducting gap)

ncp = 4e6; % n_{cp} for aluminum is 4e-6 \micro m^-3
                     % C. Wang et al. Nature Comm. 5, 5836 (2014)

Tph = 0.051; % K
tspan = [-500, 0]; % in units of \tau_0

data = load('TrapNoTrap.mat');

V_tr = data.Trap(:, 5) / delta;
% P_tr = data.Trap(:, 6);
nqp_tr = data.Trap(:, 8) - min(data.Trap(:, 8));

V_no_tr = data.NoTrap(:, 5) / delta;
% P_no_tr = data.NoTrap(:, 6);
nqp_no_tr = data.NoTrap(:, 8) - min(data.NoTrap(:, 8));

nqp_direct_no_tr = NaN(size(V_no_tr));
nqp_phonon_no_tr = NaN(size(V_no_tr));
nqp_direct_tr = NaN(size(V_tr));
nqp_phonon_tr = NaN(size(V_tr));
for k = 1:length(V_no_tr)
    if V_no_tr(k) > 1
        [~, ~, ~, ~, nqp] = directInjection0DModel(Tph, tspan,...
            V_no_tr(k), r_direct_no_tr, c_direct_no_tr);
        nqp_direct_no_tr(k) = max(nqp);
    else
        nqp_direct_no_tr(k) = 0;
    end
    if V_no_tr(k) > 2
        [~, ~, ~, ~, nqp] = phononMediatedQuasi0DModel(Tph, tspan,...
            V_no_tr(k), r_phonon_no_tr, c_phonon_no_tr);
        nqp_phonon_no_tr(k) = max(nqp);
    else
        nqp_phonon_no_tr(k) = 0;
    end
    k
end
for k = 1:length(V_tr)
    if V_tr(k) > 1
        [~, ~, ~, ~, nqp] = directInjection0DModel(Tph, tspan,...
            V_tr(k), r_direct_tr, c_direct_tr);
        nqp_direct_tr(k) = max(nqp);
    else
        nqp_direct_tr(k) = 0;
    end
    if V_tr(k) > 2
        [~, ~, ~, ~, nqp] = phononMediatedQuasi0DModel(Tph, tspan,...
            V_tr(k), r_phonon_tr, c_phonon_tr);
        nqp_phonon_tr(k) = max(nqp);
    else
        nqp_phonon_tr(k) = 0;
    end
    k
end

nqp_direct_no_tr = ncp * nqp_direct_no_tr;
nqp_direct_tr = ncp * nqp_direct_tr;
nqp_phonon_no_tr = ncp * nqp_phonon_no_tr;
nqp_phonon_tr = ncp * nqp_phonon_tr;

scrsz = get(0, 'ScreenSize');
figure('Position', [.1 .1 1.5 .8] * scrsz(4));
subplot(1, 2, 1)
[V_no_tr, nqp_no_tr, nqp_direct_no_tr, nqp_phonon_no_tr]
loglog(V_no_tr, nqp_no_tr, 'bo', V_no_tr, nqp_direct_no_tr, 'm*',...
    V_no_tr, nqp_phonon_no_tr, 'g*', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'no trap, experiment', ['direct injection: r = ',...
    num2str(r_direct_no_tr, '%.2e'), ', c = ', num2str(c_direct_no_tr, '%.2e')],...
    ['phonon-mediated excitation: r = ',...
    num2str(r_phonon_no_tr, '%.2e'), ', c = ', num2str(c_phonon_no_tr, '%.2e')]},...
    'Location', 'SouthEast')
title('No Trap')
axis tight
grid on
 
subplot(1, 2, 2)
[V_tr, nqp_tr, nqp_direct_tr, nqp_phonon_tr]
loglog(V_tr, nqp_tr, 'ko', V_tr, nqp_direct_tr, 'm*',...
    V_tr, nqp_phonon_tr, 'g*', 'MarkerSize', 10, 'LineWidth', 2)
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'trap, experiment', ['direct injection: r = ',...
    num2str(r_direct_tr, '%.2e'), ', c = ', num2str(c_direct_tr, '%.2e')],...
    ['phonon-mediated excitation: r = ',...
    num2str(r_phonon_tr, '%.2e'), ', c = ', num2str(c_phonon_tr, '%.2e')]},...
    'Location', 'SouthEast')
title('Trap')
axis tight
grid on

end