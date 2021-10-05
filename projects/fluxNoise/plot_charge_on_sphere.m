figure; hold on;
n = 400;
phi = 2*pi*rand(1,n);
costheta = -1 + 2*rand(1,n);
theta = acos(costheta);
xx = sin(theta).*cos(phi);
yy = sin(theta).*sin(phi);
zz = cos(theta);
[x y z] = sphere(100);
surface(x, y, z, 'FaceColor', 0.99*[1 1 1],'EdgeColor','none','FaceLighting','gouraud', 'FaceAlpha', 0.4)
plot3(xx(1:n/2), yy(1:n/2), zz(1:n/2), '.', 'MarkerSize', 7, 'Color', [255, 127, 14]/255)
plot3(xx(n/2:end), yy(n/2:end), zz(n/2:end), '.', 'MarkerSize', 7, 'Color', [44, 160, 44]/255)
arrow3d([0.25 -0.25], [0.25 -0.25], [-0.25 0.25], 0.85, 0.03);

daspect([1 1 1])
camlight
xl = xlim();
yl = ylim();
zl = zlim();
line(1.3*xl, [0,0], [0,0], 'LineWidth', 2, 'Color', 'k');
line([0,0], 1.3*yl, [0,0], 'LineWidth', 2, 'Color', 'k');
line([0,0], [0,0], 1.3*zl, 'LineWidth', 2, 'Color', 'k');
set(gca,'visible','off')
view(30,15)

w = 7.2/3; h = 7.2/3;
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [w h]);
set(gca, 'Position', [0 0 1 1])
set(gcf,'units','inches','position',[3,3,w,h]);
fig_path = '/Volumes/smb/mcdermott-group/users/ChrisWilen/FluxNoise/figs/';
% saveas(gcf, [fig_path 'dipole.png'])
print(gcf,[fig_path 'dipole.png'],'-dpng','-r1000');

function [h]=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The function plotting 3-dimensional arrow
%
% h=arrow3d(x,y,z,head_frac,radii,radii2,colr)
%
% The inputs are:
%       x,y,z =  vectors of the starting point and the ending point of the
%           arrow, e.g.:  x=[x_start, x_end]; y=[y_start, y_end];z=[z_start,z_end];
%       head_frac = fraction of the arrow length where the head should  start
%       radii = radius of the arrow
%       radii2 = radius of the arrow head (defult = radii*2)
%       colr =   color of the arrow, can be string of the color name, or RGB vector  (default='blue')
%
% The output is the handle of the surfaceplot graphics object.
% The settings of the plot can changed using: set(h, 'PropertyName', PropertyValue)
%
% example #1:
%        arrow3d([0 0],[0 0],[0 6],.5,3,4,[1 0 .5]);
% example #2:
%        arrow3d([2 0],[5 0],[0 -6],.2,3,5,'r');
% example #3:
%        h = arrow3d([1 0],[0 1],[-2 3],.8,3);
%        set(h,'facecolor',[1 0 0])
% 
% Written by Moshe Lindner , Bar-Ilan University, Israel.
% July 2010 (C)
if nargin==5
    radii2=radii*2;
    colr='blue';
elseif nargin==6
    colr='blue';
end
if size(x,1)==2
    x=x';
    y=y';
    z=z';
end
x(3)=x(2);
x(2)=x(1)+head_frac*(x(3)-x(1));
y(3)=y(2);
y(2)=y(1)+head_frac*(y(3)-y(1));
z(3)=z(2);
z(2)=z(1)+head_frac*(z(3)-z(1));
r=[x(1:2)',y(1:2)',z(1:2)'];
N=50;
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;
normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;
X1=[];Y1=[];Z1=[];
j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii*cos(theta).*(P1-Pc) + radii*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(2:3,j)=R1(:,1);
    Y1(2:3,j)=R1(:,2);
    Z1(2:3,j)=R1(:,3);
    j=j+1;
end
r=[x(2:3)',y(2:3)',z(2:3)'];
dr=diff(r);
dr(end+1,:)=dr(end,:);
origin_shift=(ones(size(r))*(1+max(abs(r(:))))+[dr(:,1) 2*dr(:,2) -dr(:,3)]);
r=r+origin_shift;
normdr=(sqrt((dr(:,1).^2)+(dr(:,2).^2)+(dr(:,3).^2)));
normdr=[normdr,normdr,normdr];
dr=dr./normdr;
Pc=r;
n1=cross(dr,Pc);
normn1=(sqrt((n1(:,1).^2)+(n1(:,2).^2)+(n1(:,3).^2)));
normn1=[normn1,normn1,normn1];
n1=n1./normn1;
P1=n1+Pc;
j=1;
for theta=([0:N])*2*pi./(N);
    R1=Pc+radii2*cos(theta).*(P1-Pc) + radii2*sin(theta).*cross(dr,(P1-Pc)) -origin_shift;
    X1(4:5,j)=R1(:,1);
    Y1(4:5,j)=R1(:,2);
    Z1(4:5,j)=R1(:,3);
    j=j+1;
end
X1(1,:)=X1(1,:)*0 + x(1);
Y1(1,:)=Y1(1,:)*0 + y(1);
Z1(1,:)=Z1(1,:)*0 + z(1);
X1(5,:)=X1(5,:)*0 + x(3);
Y1(5,:)=Y1(5,:)*0 + y(3);
Z1(5,:)=Z1(5,:)*0 + z(3);
h=surf(X1,Y1,Z1,'facecolor',colr,'edgecolor','none');
camlight
lighting phong

end