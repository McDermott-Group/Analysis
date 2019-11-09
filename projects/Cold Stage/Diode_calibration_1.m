% Author:   Ameya Riswadkar    
% Date:     10/7/2019
data1=loadMeasurementData('Z:\mcdermott-group\data\Cold Stage\Alexs junctions\Diode Calibration\cos2241bcz_5Mohm_1kohm_100Hz_9.hdf5');
Vdiode=data1.Current*5e6;
Idiode=data1.Voltage/1e3;
%Find offset in voltage where current is 0 for Diode
c=0;
[~,m]=min(abs(Vdiode));
c=Idiode(m);
Idiode=Idiode-c;
%plot(Vdiode,Idiode)
data2=loadMeasurementData;
data2.Current=data2.Current*5e6;
%Find offset in voltage where current is 0 for junction
[~,m]=min(abs(data2.Current));
c=data2.Voltage(m);
data2.Voltage=data2.Voltage-c;
data2.Current = interp1(Vdiode,Idiode,data2.Current);
figure
scatter(data2.Voltage, data2.Current,'filled')
xlabel('Voltage')
ylabel('Current')
grid on
hold on
index=find(data2.Voltage>=4e-4);
P = polyfit(data2.Voltage(index),data2.Current(index),1);
yfit = P(1)*data2.Voltage+P(2);
fprintf('Resistance=%f \n',1/P(1))
hold on;
plot(data2.Voltage,yfit,'r');
ax=gca;
set(gca,'FontSize',18);
xl=xlim;
yl=ylim;
x_coord=xl(1)+(xl(2)-xl(1))*0.13;
y_coord=yl(1)+(yl(2)-yl(1))*0.85;
txt=text(x_coord,y_coord,['Resistance=',num2str(1/P(1)),'\Omega'],'FontSize',18);
hold off

figure
%For Jn2
%data2.Current=data2.Current+10^(-7.289)/2;

data2.Current=data2.Current+10^(-7.3)/2;
plot(data2.Voltage, log10(abs(data2.Current)))
ax=gca;
set(gca,'FontSize',18);
xlabel('Voltage')
ylabel('Log(|Current|)')
grid on