function [r,z,charge_mat,int_val] = simulate_InducedChargeRings()
%All units in microns

%Generate the coordinate matrix for the charge placement
z_step = 1;
z_max = 375;%375 um
r_step = 3;
r_max = 550;
z = 0:z_step:z_max;
r = 1:r_step:r_max;
num_r = length(r);
num_z = length(z);
charge_mat = zeros(num_z,num_r);

%Parameters for the device rings
r_inner_step = 1;
r_inner = 108;
r_outer = 500;
r_outer_step = 1;
r_outer_max = 550;
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
%Alias the data
charge_mat = abs(mod(charge_mat + 0.5,1) - 0.5);%NOTE: Abs is not necessary

figure();imagesc(r,z,charge_mat)
title(strcat(['r_{in}: ',num2str(r_inner),', r_{out}: ',num2str(r_outer)]))
map = [zeros(10,3); parula(40)];
colormap(map)
line([0 r_inner],[0 0],'Color','red','LineWidth',4)
line([r_outer 550],[0 0],'Color','red','LineWidth',4)

[h]=makeHistogram(r,z,charge_mat);
figure(10);hold on;histogram(h,75)


axis([0.1 0.5 0 10^4])
set(gca, 'YScale', 'log')
int_val = sum(sum(h>0.1));

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