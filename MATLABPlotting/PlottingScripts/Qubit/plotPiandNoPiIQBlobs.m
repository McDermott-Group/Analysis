function [data] = plotPiandNoPiIQBlobs(gateI, gateX)
%plotPiandNoPiIQBlobs(DATA_I, DATA_X) plot IQ blobs from many trials
%
% Functions to plot two identical experiments where the difference is that
%   one of the datasets is from applying an idle gate and the other is from
%   applying a pi gate.
%
% Plots the dots and also returns stacked data in a struct.

%% Initialize and stack the data
data.Is.gateI = [];
data.Qs.gateI = [];
data.Is.gateX = [];
data.Qs.gateX = [];

for n = 1:length(gateI.Trial)
   data.Is.gateI = [data.Is.gateI gateI.Is(n,:)];
   data.Qs.gateI = [data.Qs.gateI gateI.Qs(n,:)];
   data.Is.gateX = [data.Is.gateX gateX.Is(n,:)];
   data.Qs.gateX = [data.Qs.gateX gateX.Qs(n,:)];
end

data.Is.all = [data.Is.gateI(1,:) data.Is.gateX(1,:)]';
data.Qs.all = [data.Qs.gateI(1,:) data.Qs.gateX(1,:)]';

data.allIQs = [data.Is.all data.Qs.all];

%% Calculate mean IQ values via real data and via kmeans clustering
% mean values
data.I.gateI = mean(data.Is.gateI);
data.Q.gateI = mean(data.Qs.gateI);

data.I.gateX = mean(data.Is.gateX);
data.Q.gateX = mean(data.Qs.gateX);

% kmeans clustering
opts = statset('Display','final',...
               'UseParallel', 1,...
               'MaxIter', 10000);
[data.idx, data.C] = kmeans(data.allIQs, 2, 'Replicates', 1,...
                                            'Options', opts);

%% Do the plotting
% mean values
figure;
subplot(1,2,1)
plot(data.Is.gateI, data.Qs.gateI,'b.',...
    'MarkerSize',2,'DisplayName','0')
hold on
plot(data.Is.gateX, data.Qs.gateX,'r.',...
    'MarkerSize',2,'DisplayName','1')
xlabel('I (V)')
ylabel('Q (V)')
title({'Raw Quadrature Space Plot for ', ...
    [num2str(length(data.Is.gateI)) ' Shots']})
axis square
grid on
set(gca,'FontSize',14)
plot([data.I.gateI data.I.gateX],[data.Q.gateI data.Q.gateX],'g+',...
    'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');
legend('show')
hold off

% kmeans clustering
if data.idx(1) == 1
    blobAColor = 'b.';
    blobBColor = 'r.';
else
    blobAColor = 'r.';
    blobBColor = 'b.';
end
subplot(1,2,2)
plot(data.allIQs(data.idx==1,1),data.allIQs(data.idx==1,2),blobAColor,...
    'MarkerSize',2,'DisplayName','Blob A')
hold on
plot(data.allIQs(data.idx==2,1),data.allIQs(data.idx==2,2),blobBColor,...
    'MarkerSize',2,'DisplayName','Blob B')
plot(data.C(:,1),data.C(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3,'DisplayName','Fitted Centers')
xlabel('I (V)')
ylabel('Q (V)')
title({'Clustered Quadrature Space Plot for ',...
    [num2str(length(data.Is.gateI)) ' Shots']})
axis square
grid on
set(gca,'FontSize',14)
plot([data.I.gateI data.I.gateX],[data.Q.gateI data.Q.gateX],'g+',...
    'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');

legend('show')
end