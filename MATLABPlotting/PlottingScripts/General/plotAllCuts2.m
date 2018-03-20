function plotAllCuts2
% Function to take an existing data plot that's already open in MATLAB (and
% in focus) and plots all cuts along the Y axis on a common plot instead of
% using the 2D false color plot.
%
% Note: if you aren't sure what plot is in focus and want to be sure, type 
%
%      figure(n)
%
% where n is the figure number assigned to the figure you wish to cut.

h = gcf;
axesObjs = get(h, 'Children');
dataObjs = get(axesObjs, 'Children');
temp = dataObjs{2};
cD = get(temp, 'CData');
yD = get(temp, 'YData');
xD = get(temp, 'XData');
figure();
hold on
for n=1:length(cD(:,1))
    plot(xD, cD(n,:), 'DisplayName', num2str(yD(n)))
end
legend(gca,'show')
title(h.CurrentAxes.Title.String);
xlabel(h.CurrentAxes.XLabel.String);
ylabel([h.CurrentAxes.Title.String{1}]);
set(gca,'FontSize',14);
end