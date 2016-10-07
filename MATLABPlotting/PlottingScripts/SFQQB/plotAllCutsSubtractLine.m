function [  ] = plotAllCutsSubtractLine(m,b)
% Function to take an existing data plot that's already open in MATLAB (and
% in focus) and plots all cuts along the X axis on a common plot instead of
% using the 2D false color plot.
%
% n.b. if you aren't sure what plot is in focus and want to be sure, type 
%
%      figure(n)
%
% where n is the figure number assigned to the figure you wish to cut.

h = gcf;
axesObjs = get(h,'Children'); dataObjs = get(axesObjs,'Children');
temp = dataObjs{2};
cD = get(temp,'CData');
yD = get(temp,'YData');
xD = get(temp,'XData');

for j = 1:length(xD)
   for k = 1:length(yD)
      cD(k,j) =  cD(k,j) - m * xD(j) - b;
   end
end

figure();
hold on
for n=1:length(cD(1,:))
    plot(yD,cD(:,n),'DisplayName',num2str(xD(n)))
    
end
legend(gca,'show')
title(h.CurrentAxes.Title.String);
xlabel(h.CurrentAxes.YLabel.String);
ylabel([h.CurrentAxes.Title.String{1}]);
set(gca,'FontSize',14);
end