h1 = errorbar(IdleTime/1000,GateFidelity_DissOFF,GateFidelity_DissOFF_err);
hold on
h2 = errorbar(IdleTimeOn/1000,GateFidelity_DissON,GateFidelity_DissON_err);
% t1 = plot(IdleTime/1000, T1limit);
% title({'Idle Gate'},'interpreter','latex','FontSize',16)

xlabel({'Idle Gate Time (usec)'},'interpreter','latex','FontSize',16)

ylabel({'Fidelity'},'interpreter','latex','FontSize',16)

h1.LineStyle = 'none'; 
h2.LineStyle = 'none'; 
h1.Marker = '.';
h1.Color = [0.635 0.078 0.184];
h2.Color = [0 0.447 0.741];
h2.Marker = '.'; 
h1.MarkerSize = 20;
h2.MarkerSize = 20;
% t1.LineWidth = 1;
% t1.Color = 'b';
% leg = legend( [h1,h2,t1], 'Dissipator OFF', 'Dissipator ON','T1 limit', 'Location', 'NorthEast' );
hold off
grid on
axis tight