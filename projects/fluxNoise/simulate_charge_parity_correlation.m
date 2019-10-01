N = 1000;
nFiles = 100;
dt = 100e-6;
T_qp = 1e-3;
T1 = 70e-6;
M_time = 3e-6;
correlation = 1/10;

r = rand(N, nFiles);
dP = r < 1/T_qp*dt;
P = mod(cumsum(dP,1),2);

r = rand(N, nFiles);
M_cor = r > dP*correlation;

r = rand(N, nFiles);
M_decay = r > exp(-M_time/T1);

M = M_cor & M_decay;

tot = zeros(N,1);
t = (1:N)*dt;
for i = 1:nFiles
    
    for j = 1:N
        tot(j) = tot(j) + sum( P(1:end-j,i)==P(j+1:end,i) )/(N-j)/nFiles;
    end
end

figure; plot(t,tot-0.5); grid on;



sum11ee = zeros(N,1);
sum11eo = zeros(N,1);
sum10ee = zeros(N,1);
sum10eo = zeros(N,1);
sum01ee = zeros(N,1);
sum01eo = zeros(N,1);
sum00ee = zeros(N,1);
sum00eo = zeros(N,1);
t = (1:N)*100e-6;
for i = 1:nFiles
    
    for j = 1:N
        dP = P(1:end-j,i)==P(j+1:end,i);
        Ma = M(1:end-j,i);  Mb = M(j+1:end,i);
        M0a = (Ma == 0);  M0b = (Mb == 0);
        M1a = (Ma == 1);  M1b = (Mb == 1);
        sum11ee(j) = sum11ee(j) + sum((M1a==1) + (M1b==1) + (dP==1) == 3)/j/nFiles;
        sum11eo(j) = sum11eo(j) + sum((M1a==1) + (M1b==1) + (dP==0) == 3)/j/nFiles;
        sum10ee(j) = sum10ee(j) + sum((M1a==1) + (M1b==0) + (dP==1) == 3)/j/nFiles;
        sum10eo(j) = sum10eo(j) + sum((M1a==1) + (M1b==0) + (dP==0) == 3)/j/nFiles;
        sum01ee(j) = sum01ee(j) + sum((M1a==0) + (M1b==1) + (dP==1) == 3)/j/nFiles;
        sum01eo(j) = sum01eo(j) + sum((M1a==0) + (M1b==1) + (dP==0) == 3)/j/nFiles;
        sum00ee(j) = sum00ee(j) + sum((M1a==0) + (M1b==0) + (dP==1) == 3)/j/nFiles;
        sum00eo(j) = sum00eo(j) + sum((M1a==0) + (M1b==0) + (dP==0) == 3)/j/nFiles;
    end
end

figure; plot(t,[sum11ee,sum11eo,sum10ee,sum10eo,sum01ee,sum01eo,sum00ee,sum00eo]); grid on;
legend('R_{11}^{ee}','R_{11}^{eo}','R_{10}^{ee}','R_{10}^{eo}','R_{01}^{ee}','R_{01}^{eo}','R_{00}^{ee}','R_{00}^{eo}');