function [r,z,charge_mat,int_val] = simulate_InducedChargeRings()
%All units in microns

%%%%%%%% for figure 1 %%%%%%%%%
% w = 8.9/2; h = 0.65*2.54;

% %Generate the coordinate matrix for the charge placement
% z_step = 1;
% z_max = 375;%375 um
% r_step = 1;
% r_max = 200;
% z = 0:z_step:z_max;
% r = 1:r_step:r_max;
% num_r = length(r);
% num_z = length(z);
% charge_mat = zeros(num_z,num_r);
% 
% %Parameters for the device rings
% r_inner_step = 1;
% r_inner = 70;
% r_outer = 90.5;
% r_outer_step = 1;
% r_outer_max = 550;
% z_shift = 0.1;%To avoid Infs

%%%%%%%% for expanded range %%%%%%%%%
w = 8.9; h = 2*0.65*2.54;

%Generate the coordinate matrix for the charge placement
z_step = 1;
z_max = 375;%375 um
r_step = 3;
r_max = 1000;
z = 0:z_step:z_max;
r = 1:r_step:r_max;
num_r = length(r);
num_z = length(z);
charge_mat = zeros(num_z,num_r);

%Parameters for the device rings
r_inner_step = 1.5;
r_inner = 70;
r_outer = 90.5;
r_outer_step = 4;
r_outer_max = 1200;
z_shift = 0.1;%To avoid Infs

%Initialize vectors
r_coord_vec = [r_inner_step:r_inner_step:r_inner, r_outer:r_outer_step:r_outer_max,0]';
z_coord_vec = zeros(size(r_coord_vec));
r_coord_num = length(r_coord_vec);
r_inner_num = length(r_inner_step:r_inner_step:r_inner);
v_vec = zeros(size(r_coord_vec));
v_vec(end) = 1;

%Precompute C_mat_inv(1:n-1, 1:n-1)
Cinv_mat = zeros(r_coord_num,r_coord_num);
for i = 1:(r_coord_num-1)
    Cinv_mat(:,i) = potential(r_coord_vec(i),z_coord_vec(i),r_coord_vec,z_coord_vec,z_shift);
end

for zi = 1:num_z
    z_coord_vec(end) = z(zi);
    for ri = 1:num_r
        r_coord_vec(end) = r(ri);

        % Cinv_mat = zeros(r_coord_num,r_coord_num);
        % for i = 1:r_coord_num
        %     Cinv_mat(:,i) = potential(r_coord_vec(i),z_coord_vec(i),r_coord_vec,z_coord_vec,z_shift);
        % end
        last_set = potential(r_coord_vec(end),z_coord_vec(end),r_coord_vec,z_coord_vec,z_shift);
        Cinv_mat(:,end) = last_set;
        Cinv_mat(end,:) = last_set';
        
        q_vec = Cinv_mat\v_vec;% inv(Cinv_mat)*v_vec
        charge_inner = sum(q_vec(1:r_inner_num));
        % charge_outer = sum(q_vec((r_inner_num+1):(end-1)));
        charge = q_vec(end);
        
        charge_mat(zi,ri) = charge_inner/charge;
    end
end
charge_mat(isnan(charge_mat)) = -1;
%Alias the data
% charge_mat = abs(noiselib.alias(charge_mat, 0.5));%NOTE: Abs is not necessary
% Plot Gradient
[gz, gr] = gradient(charge_mat, z_step*1e-6, r_step*1e-6);
alpha_grad = sqrt(gz.^2 + gr.^2);

plot_and_save(-charge_mat, r, z, r_inner, r_outer, w, h, 'alpha.pdf')
plot_and_save(alpha_grad, r, z, r_inner, r_outer, w, h, 'alpha_gradient.pdf')

[h]=makeHistogram(r,z,charge_mat);
figure(10);hold on;histogram(h,75)

axis([0.1 0.5 0 10^4])
set(gca, 'YScale', 'log')
int_val = sum(sum(h>0.1));

end


function [] = plot_and_save(A, r, z, r1, r2, w, h, file)

figure(); imagesc([flip(-r) r],z,[flip(A,2) A])
%title(strcat(['r_{in}: ',num2str(r1),', r_{out}: ',num2str(r2)]))
map = [zeros(0,3); parula(50)];
colormap(map)
rectangle('Position', [-r1   -20 2*r1    20], 'FaceColor', 'red', 'LineStyle', 'none');
rectangle('Position', [ r2   -20 1500-r2 20], 'FaceColor', 'red', 'LineStyle', 'none');
rectangle('Position', [-1500 -20 1500-r2 20], 'FaceColor', 'red', 'LineStyle', 'none');
ylim([-20,inf]);

axis off;
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [w h]);
set(gca, 'Position', [0 0 1 1])
cb = colorbar('south');%,'color','white');
pos = cb.Position;
cb.Position = [0.5-0.035,pos(2),pos(3)/2-0.035,2*pos(4)];
set(gcf,'units','centimeters','position',[3,3,w,h]);
set(cb,'FontSize',10);
set(cb,'FontName','Arial');
set(gca,'ColorScale','log')
l = floor(log10(min(A(:))));
u = ceil(log10(max(A(:))));
caxis(10.^[ l, u ]);
set(cb,'YTick',10.^flip(u:-floor((u-l)/2):l))
fig_path = '/Volumes/smb/mcdermott-group/users/ChrisWilen/FluxNoise/figs/';
%saveas(gcf, [fig_path file])

end

function[v] = potential(r,zr,p,zp,z_shift)
%Calculate the potential at a point (p,0,zp) due to a ring of charge Q=1
%centered at (0,0,zr) with radius r.
z = zr - zp + z_shift;
rd = r - p;
m = -4*r.*p./(z.^2+rd.^2);
ellk = 1./sqrt(1-m).*ellipke(-m./(1-m));% Imaginary-modulus transformation
v = 1/(2*pi) * 4./sqrt(z.^2+rd.^2).*ellk;
end

function[h] = makeHistogram(r,z,charge_mat)
h = zeros(length(z),sum(r));
ind = 1;
for i = 1:length(r)
    for j = 1:r(i)
        ind = ind + 1;
        h(:,ind) = charge_mat(:,i);
    end
end
end