function mixedFit2TrapNoTrap
%mixedFit2TrapNoTrap Fitting to TrapNoTrap dataset.

r_direct_no_tr = 5e-6; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon_no_tr = 1e-10; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_no_tr = 0; % trapping rate in units of 1 / \tau_0

r_direct_tr = 5e-6; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
r_phonon_tr = 1e-10; % in units of 1 / \tau_0 %(assuming n_{qp} in units of n_{cp})
c_tr = .1; % trapping rate in units of 1 / \tau_0

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

nqp_sim_no_tr = NaN(size(V_no_tr));
nqp_sim_tr = NaN(size(V_tr));
for k = 1:length(V_no_tr)
    if V_no_tr(k) > 1
        [~, ~, ~, ~, nqp] = mixedInjectionQuasi0DModel(Tph, tspan,...
            V_no_tr(k), r_direct_no_tr, r_phonon_no_tr, c_no_tr);
        nqp_sim_no_tr(k) = max(nqp);
    else
        nqp_sim_no_tr(k) = 0;
    end
    k
end
for k = 1:length(V_tr)
    if V_tr(k) > 1
        [~, ~, ~, ~, nqp] = mixedInjectionQuasi0DModel(Tph, tspan,...
            V_tr(k), r_direct_tr, r_phonon_tr, c_tr);
        nqp_sim_tr(k) = max(nqp);
    else
        nqp_sim_tr(k) = 0;
    end
    k
end

nqp_sim_no_tr = ncp * nqp_sim_no_tr;
nqp_sim_tr = ncp * nqp_sim_tr;

scrsz = get(0, 'ScreenSize');
figure('Position', [.1 .1 1.5 .8] * scrsz(4));
subplot(1, 2, 1)
loglog(V_no_tr, nqp_no_tr, 'bo', V_no_tr, nqp_sim_no_tr, 'm*',...
    'MarkerSize', 10, 'LineWidth', 2)
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'no trap, experiment', ['simulation: r_{qp} = ',...
    num2str(r_direct_no_tr, '%.2e'), ', r_{ph} = ',...
    num2str(r_phonon_no_tr, '%.2e'), ', c = ',...
    num2str(c_no_tr, '%.2e')]},...
    'Location', 'SouthEast')
title('No Trap')
axis tight
grid on
 
subplot(1, 2, 2)
loglog(V_tr, nqp_tr, 'ko', V_tr, nqp_sim_tr, 'm*',...
    'MarkerSize', 10, 'LineWidth', 2)
xlabel('Injection Energy (\Delta)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu m^{-3})', 'FontSize', 14)
legend({'experiment', ['simulation: r_{qp} = ',...
    num2str(r_direct_tr, '%.2e'), ', r_{ph} = ',...
    num2str(r_phonon_tr, '%.2e'), ', c = ',...
    num2str(c_tr, '%.2e')]},...
    'Location', 'SouthEast')
title('Trap')
axis tight
grid on

end