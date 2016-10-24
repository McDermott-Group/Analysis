function [data] = plotPiandNoPiIQBlobs(gateI, gateX)
%plotPiandNoPiIQBlobs(DATA_I, DATA_X) plot IQ blobs from many trials
%
% Functions to plot two identical experiments where the difference is that
%   one of the datasets is from applying an idle gate and the other is from
%   applying a pi gate.
%
% Plots the dots and also returns stacked data in a struct.

%% Initialize and stack the data
data.gateI.Is = [];
data.gateI.Qs = [];
data.gateX.Is = [];
data.gateX.Qs = [];

for n = 1:length(gateI.Trial)
   data.gateI.Is = [data.gateI.Is gateI.Is(n,:)];
   data.gateI.Qs = [data.gateI.Qs gateI.Qs(n,:)];
   data.gateX.Is = [data.gateX.Is gateX.Is(n,:)];
   data.gateX.Qs = [data.gateX.Qs gateX.Qs(n,:)];
end

data.gateI.IQs = [data.gateI.Is' data.gateI.Qs'];
data.gateX.IQs = [data.gateX.Is' data.gateX.Qs'];

data.all.Is = [data.gateI.Is(1,:) data.gateX.Is(1,:)]';
data.all.Qs = [data.gateI.Qs(1,:) data.gateX.Qs(1,:)]';
data.all.IQs = [data.all.Is data.all.Qs];

%% Calculate mean IQ values via real data and via kmeans clustering
% mean values
data.gateI.I = mean(data.gateI.Is);
data.gateI.Q = mean(data.gateI.Qs);

data.gateX.I = mean(data.gateX.Is);
data.gateX.Q = mean(data.gateX.Qs);

% kmeans clustering
opts = statset('Display','final',...
               'UseParallel', 1,...
               'MaxIter', 10000);
[data.gateI.idx, data.gateI.C] = kmeans(data.all.IQs, 2,...
                                            'Replicates', 1,...
                                            'Options', opts);          

[data.all.idx, data.all.C] = kmeans(data.all.IQs, 2, 'Replicates', 1,...
                                            'Options', opts);

%% Do the plotting
% mean values
figure;
subplot(1,2,1)
plot(data.gateI.Is, data.gateI.Qs,'b.',...
    'MarkerSize',2,'DisplayName','0')
hold on
plot(data.gateX.Is, data.gateX.Qs,'r.',...
    'MarkerSize',2,'DisplayName','1')
xlabel('I (V)')
ylabel('Q (V)')
title({'Raw Quadrature Space Plot for ', ...
    [num2str(length(data.gateI.Is)) ' Shots']})
axis square
grid on
set(gca,'FontSize',14)
plot([data.gateI.I data.gateX.I],[data.gateI.Q data.gateX.Q],'g+',...
    'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');
legend('show')
hold off

% kmeans clustering
if data.all.idx(1) == 1
    blobAColor = 'b.';
    blobBColor = 'r.';
else
    blobAColor = 'r.';
    blobBColor = 'b.';
end
subplot(1,2,2)
plot(data.all.IQs(data.all.idx==1,1),data.all.IQs(data.all.idx==1,2),blobAColor,...
    'MarkerSize',2,'DisplayName','Blob A')
hold on
plot(data.all.IQs(data.all.idx==2,1),data.all.IQs(data.all.idx==2,2),blobBColor,...
    'MarkerSize',2,'DisplayName','Blob B')
plot(data.all.C(:,1),data.all.C(:,2),'kx',...
    'MarkerSize',15,'LineWidth',3,'DisplayName','Fitted Centers')
xlabel('I (V)')
ylabel('Q (V)')
title({'Clustered Quadrature Space Plot for ',...
    [num2str(length(data.gateI.Is)) ' Shots']})
axis square
grid on
set(gca,'FontSize',14)
plot([data.gateI.I data.gateX.I],[data.gateI.Q data.gateX.Q],'g+',...
    'MarkerSize',10,'LineWidth',3,'DisplayName','Actual Blob Centers');

legend('show')
end