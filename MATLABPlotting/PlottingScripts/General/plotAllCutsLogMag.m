function [ h ] = plotAllCutsLogMag()
h = gcf;
axesObjs = get(h,'Children'); dataObjs = get(axesObjs,'Children');
temp = dataObjs{2};
cDlogMag = 20*log10(get(temp,'CData'));
yD = get(temp,'YData');
xD = get(temp,'XData');
figure();
hold on
for n=1:length(cDlogMag(1,:))
    plot(yD,cDlogMag(:,n),'DisplayName',num2str(xD(n)))
    
end
legend(gca,'show')
title(h.CurrentAxes.Title.String);
xlabel(h.CurrentAxes.YLabel.String);
ylabel([h.CurrentAxes.Title.String{1}, 'LogMag']);
set(gca,'FontSize',14);
end