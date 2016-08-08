function steadyStateModelFit2NqpWorkSpace
% steadyStateModelFit2NqpWorkSpace Fitting to the `NqpWorkSpace` dataset
% using the two-region steady-state quasi-0D model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near, Cuts, No Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 6.508e-05; r_phonon = 1.467e+00; c = 2.335e-02; vol = 1.306e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.CutNear(:, 1) / delta;
P = data.CutNear(:, 2);
nqp = data.CutNear(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'b^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceCutNearV.pdf')
print(h, 'SimNqpWorkSpaceCutNearV.fig')

h = figure;
hold on
plot(P, nqp, 'b^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceCutNearP.pdf')
print(h, 'SimNqpWorkSpaceCutNearP.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near, No Cuts, No Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 5.519e-05; r_phonon = 1.472e+00; c = 2.397e-02; vol = 1.572e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.NoCutNear(:, 1) / delta;
P = data.NoCutNear(:, 2);
nqp = data.NoCutNear(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, No Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceNoCutNearV.pdf')
print(h, 'SimNqpWorkSpaceNoCutNearV.fig')

h = figure;
hold on
plot(P, nqp, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, No Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceNoCutNearP.pdf')
print(h, 'SimNqpWorkSpaceNoCutNearP.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Near, Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 7.397e-06; r_phonon = 1.786e-01; c = 6.825e-02; vol = 2.034e+05;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.TrapNear(:, 1) / delta;
P = data.TrapNear(:, 2);
nqp = data.TrapNear(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceTrapNearV.pdf')
print(h, 'SimNqpWorkSpaceTrapNearV.fig')

h = figure;
hold on
plot(P, nqp, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Near Resonator, Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceTrapNearP.pdf')
print(h, 'SimNqpWorkSpaceTrapNearP.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far, Cuts, No Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 2.771e-05; r_phonon = 9.358e-01; c = 2.854e-02; vol = 3.293e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.CutFar(:, 1) / delta;
P = data.CutFar(:, 2);
nqp = data.CutFar(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'b^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceCutFarV.pdf')
print(h, 'SimNqpWorkSpaceCutFarV.fig')

h = figure;
hold on
plot(P, nqp, 'b^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceCutFarP.pdf')
print(h, 'SimNqpWorkSpaceCutFarP.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far, No Cuts, No Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 2.411e-05; r_phonon = 9.481e-01; c = 2.804e-02; vol = 3.787e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.NoCutFar(:, 1) / delta;
P = data.NoCutFar(:, 2);
nqp = data.NoCutFar(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, No Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceNoCutFarV.pdf')
print(h, 'SimNqpWorkSpaceNoCutFarV.fig')

h = figure;
hold on
plot(P, nqp, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, No Cuts, No Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceNoCutFarP.pdf')
print(h, 'SimNqpWorkSpaceNoCutFarP.fig')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far, Traps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 3.135e-05; r_phonon = 1.020e+00; c = 8.442e-02; vol = 4.017e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

% Data.
delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.TrapFar(:, 1) / delta;
P = data.TrapFar(:, 2);
nqp = data.TrapFar(:, 4);

V_sim = V;
nqp_sim = NaN(size(V_sim));
P_sim = NaN(size(V_sim));
parfor k = 1:length(V_sim)
    if V_sim(k) > 1
        [~, ~, ~, ~, nqp_sim_t, ~, P_sim(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim(k), r_direct, r_phonon, c, c, vol, N);
        nqp_sim(k) = max(nqp_sim_t);
    else
        nqp_sim(k) = 0;
    end
    fprintf('*')
end

h = figure;
hold on
plot(V, nqp, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(V_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceTrapFarV.pdf')
print(h, 'SimNqpWorkSpaceTrapFarV.fig')

h = figure;
hold on
plot(P, nqp, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2)
plot(P_sim, nqp_sim, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'm')
hold off
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment', 'simulation'},...
        'Location', 'SouthEast')
title({'Far Resonator, Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'xscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNqpWorkSpaceTrapFarP.pdf')
print(h, 'SimNqpWorkSpaceTrapFarP.fig')

end