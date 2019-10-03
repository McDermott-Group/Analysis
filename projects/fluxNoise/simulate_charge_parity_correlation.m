N = 1000;
nFiles = 100;
dt = 100e-6;        % time between runs
T_qp = 1e-3;        % QP poisoning time
T1 = 100e-6;        % T1 of qubit
M_time = 3e-6;      % time before measuring qubit
correlation = 9/10; % what portion of tunneling events causes a qubit state change

% calc when there is a parity change (dP) and the parity for each meas (P)
r = rand(N, nFiles);
dP = r < 1/T_qp*dt;
P = mod(cumsum(dP,1),2);

% for each parity change in dP, does it cause a state change (M_cor)?
r = rand(N, nFiles);
M_cor = dP.*(r < dP*correlation);

% calculate a measurement based on normal decay but not parity change
r = rand(N, nFiles);
M_decay = r < exp(-M_time/T1);

% combine state change due to parity change and normal decay
M = mod(M_decay + M_cor, 2);

t = (1:N)*dt;

% Plot autocorrelation of Parity P
tot = zeros(N,1);
for i = 1:nFiles
    for j = 1:N
        tot(j) = tot(j) + sum( P(1:end-j,i)==P(j+1:end,i) )/(N-j)/nFiles;
    end
end
figure; plot(t,tot); grid on;


% Plot autocorrelation between charge parity and qubit state
sum11ee = zeros(N,1);
sum11eo = zeros(N,1);
sum10ee = zeros(N,1);
sum10eo = zeros(N,1);
sum01ee = zeros(N,1);
sum01eo = zeros(N,1);
sum00ee = zeros(N,1);
sum00eo = zeros(N,1);
for i = 1:nFiles
    
    for j = 1:N
        dPj = P(1:end-j,i)==P(j+1:end,i);
        Ma = M(1:end-j,i);  Mb = M(j+1:end,i);
        M0a = (Ma == 0);  M0b = (Mb == 0);
        M1a = (Ma == 1);  M1b = (Mb == 1);
%         sum11ee(j) = sum11ee(j) + sum((M1a==1) + (M1b==1) + (dPj==1) == 3)/j/nFiles;
%         sum11eo(j) = sum11eo(j) + sum((M1a==1) + (M1b==1) + (dPj==0) == 3)/j/nFiles;
%         sum10ee(j) = sum10ee(j) + sum((M1a==1) + (M1b==0) + (dPj==1) == 3)/j/nFiles;
%         sum10eo(j) = sum10eo(j) + sum((M1a==1) + (M1b==0) + (dPj==0) == 3)/j/nFiles;
%         sum01ee(j) = sum01ee(j) + sum((M1a==0) + (M1b==1) + (dPj==1) == 3)/j/nFiles;
%         sum01eo(j) = sum01eo(j) + sum((M1a==0) + (M1b==1) + (dPj==0) == 3)/j/nFiles;
%         sum00ee(j) = sum00ee(j) + sum((M1a==0) + (M1b==0) + (dPj==1) == 3)/j/nFiles;
%         sum00eo(j) = sum00eo(j) + sum((M1a==0) + (M1b==0) + (dPj==0) == 3)/j/nFiles;
        sum11ee(j) = sum11ee(j) + sum((M1a==1) + (M1b==1) + (dPj==1) == 3);
        sum11eo(j) = sum11eo(j) + sum((M1a==1) + (M1b==1) + (dPj==0) == 3);
        sum10ee(j) = sum10ee(j) + sum((M1a==1) + (M1b==0) + (dPj==1) == 3);
        sum10eo(j) = sum10eo(j) + sum((M1a==1) + (M1b==0) + (dPj==0) == 3);
        sum01ee(j) = sum01ee(j) + sum((M1a==0) + (M1b==1) + (dPj==1) == 3);
        sum01eo(j) = sum01eo(j) + sum((M1a==0) + (M1b==1) + (dPj==0) == 3);
        sum00ee(j) = sum00ee(j) + sum((M1a==0) + (M1b==0) + (dPj==1) == 3);
        sum00eo(j) = sum00eo(j) + sum((M1a==0) + (M1b==0) + (dPj==0) == 3);
    end
end

sum11 = sum11ee + sum11eo;
sum10 = sum10ee + sum10eo;
sum01 = sum01ee + sum01eo;
sum00 = sum00ee + sum00eo;

% figure; plot(t,[sum11ee,sum11eo,sum10ee,sum10eo,sum01ee,sum01eo,sum00ee,sum00eo]); grid on;
% legend('R_{11}^{ee}','R_{11}^{eo}','R_{10}^{ee}','R_{10}^{eo}','R_{01}^{ee}','R_{01}^{eo}','R_{00}^{ee}','R_{00}^{eo}');
% figure; plot(t,[sum11ee./sum11,sum11eo./sum11,sum10ee./sum10,sum10eo./sum10,sum01ee./sum01,sum01eo./sum01,sum00ee./sum00,sum00eo./sum00]); grid on;
% legend('R_{11}^{ee}','R_{11}^{eo}','R_{10}^{ee}','R_{10}^{eo}','R_{01}^{ee}','R_{01}^{eo}','R_{00}^{ee}','R_{00}^{eo}');
figure; plot(t,[sum11ee./sum11,sum10ee./sum10,sum01ee./sum01,sum00ee./sum00]); grid on;
legend('R_{11}^{ee}','R_{10}^{ee}','R_{01}^{ee}','R_{00}^{ee}');