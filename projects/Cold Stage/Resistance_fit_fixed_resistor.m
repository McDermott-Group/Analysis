% Author:   Ameya Riswadkar    
% Date:     10/7/2019
data1=loadMeasurementData;
%Find offset in voltage where current is 0 for junction
[~,m]=min(abs(data1.Current));
c=data1.Voltage(m);
data1.Voltage=data1.Voltage-c;
figure
scatter(data1.Voltage, data1.Current,'filled')
xlabel('Voltage')
ylabel('Current')
grid on
hold on
index=find(data1.Voltage>=4.25e-4);
P = polyfit(data1.Voltage(index),data1.Current(index),1);
yfit = P(1)*data1.Voltage+P(2);
fprintf('Resistance=%f \n',1/P(1))
hold on;
plot(data1.Voltage,yfit,'r');
xl=xlim;
yl=ylim;
x_coord=xl(1)+(xl(2)-xl(1))*0.13;
y_coord=yl(1)+(yl(2)-yl(1))*0.85;
txt=text(x_coord,y_coord,['Resistance=',num2str(1/P(1)),'\Omega'],'FontSize',18);
ax=gca;
set(gca,'FontSize',18);
%txt.FontSize=28;
hold off