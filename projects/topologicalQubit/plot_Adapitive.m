[file,path] = uigetfile('*.*');
if isequal(file,0)
   disp('User selected Cancel');
else
   h5disp([path file])
end
amplitude = h5read([path file], '/dependents/Occupation Data');
f_qb = h5read([path file], '/dependents/Center Frequency');
flux = h5read([path file], '/independents/Plaq Flux Bias');
f = h5read([path file], '/independents/QB Drive Frequency');
figure(101)
plot(flux(1:40), f_qb(1:40))