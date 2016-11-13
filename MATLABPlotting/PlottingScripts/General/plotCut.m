% figure(2);
h = gcf;
axesObjs = get(h,'Children'); dataObjs = get(axesObjs,'Children');
temp = dataObjs{2};
cD = get(temp,'CData');
yD = get(temp,'YData');
xD = get(temp,'XData');
figure();
plot(yD,cD(:,13))