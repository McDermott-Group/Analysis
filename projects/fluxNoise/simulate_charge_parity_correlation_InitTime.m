idle_times = (10:10:1000)*1e-6;
reps = 4000;
n_files = 2;
dt = 10e-9;         % time between runs
T_qp = 1e-3;        % QP poisoning time
T1 = 20e-6;        % T1 of qubit
M_time = 3e-6;      % time before measuring qubit
correlation = 3/10; % what portion of tunneling events causes a qubit state change

P  = zeros(n_files,length(idle_times),reps);
M1 = zeros(n_files,length(idle_times),reps);
M2 = zeros(n_files,length(idle_times),reps);
M1(:,:,1) = 1; M2(:,:,1) = 1;
for i = 1:n_files
    for j = 1:length(idle_times)
        for k = 2:reps
            [M1(i,j,k), P(i,j,k)] = evolve_state(M2(i,j,k-1), P(i,j,k-1), idle_times(j), T_qp, correlation, T1, dt);
            M2(i,j,k) = mod(M1(i,j,k)+P(i,j,k),2);
        end
    end
end

% check to make sure we see T1 decay in M and T_qp decay in P
% avgM = zeros(length(idle_times),1);
% avgP = zeros(length(idle_times),1);
% for i = 1:length(idle_times)
%     for j = 1:100
%         [m,p] = evolve_state(1, 1, idle_times(i), T_qp, 0, T1, dt);
%         avgM(i) = avgM(i) + m/100;
%         avgP(i) = avgP(i) + p/100;
%     end
% end
% plot(idle_times,avgM,idle_times,avgP)

% 
% % Plot autocorrelation of Parity P
% tot = zeros(N,1);
% for i = 1:n_files
%     for j = 1:N
%         tot(j) = tot(j) + sum( P(1:end-j,i)==P(j+1:end,i) )/(N-j)/n_files;
%     end
% end
% figure; plot(t,tot); grid on;


% Plot autocorrelation between charge parity and qubit state
P = M1 == M2;
dP = abs(P(:,:,2:end) - P(:,:,1:end-1));
Ma = M2(:,:,1:end-1);
Mb = M1(:,:,2:end);
sum11ee = sum(sum((Ma==1) + (Mb==1) + (dP==0) == 3, 3));
sum11eo = sum(sum((Ma==1) + (Mb==1) + (dP==1) == 3, 3));
sum10ee = sum(sum((Ma==1) + (Mb==0) + (dP==0) == 3, 3));
sum10eo = sum(sum((Ma==1) + (Mb==0) + (dP==1) == 3, 3));
sum01ee = sum(sum((Ma==0) + (Mb==1) + (dP==0) == 3, 3));
sum01eo = sum(sum((Ma==0) + (Mb==1) + (dP==1) == 3, 3));
sum00ee = sum(sum((Ma==0) + (Mb==0) + (dP==0) == 3, 3));
sum00eo = sum(sum((Ma==0) + (Mb==0) + (dP==1) == 3, 3));


sum11 = sum11ee + sum11eo;
sum10 = sum10ee + sum10eo;
sum01 = sum01ee + sum01eo;
sum00 = sum00ee + sum00eo;

% figure; plot(t,[sum11ee,sum11eo,sum10ee,sum10eo,sum01ee,sum01eo,sum00ee,sum00eo]); grid on;
% legend('R_{11}^{ee}','R_{11}^{eo}','R_{10}^{ee}','R_{10}^{eo}','R_{01}^{ee}','R_{01}^{eo}','R_{00}^{ee}','R_{00}^{eo}');
% figure; plot(t,[sum11ee./sum11,sum11eo./sum11,sum10ee./sum10,sum10eo./sum10,sum01ee./sum01,sum01eo./sum01,sum00ee./sum00,sum00eo./sum00]); grid on;
% legend('R_{11}^{ee}','R_{11}^{eo}','R_{10}^{ee}','R_{10}^{eo}','R_{01}^{ee}','R_{01}^{eo}','R_{00}^{ee}','R_{00}^{eo}');
figure; plot(idle_times,[(sum11ee - sum11eo)'./sum11',(sum10ee - sum10eo)'./sum10',(sum01ee - sum01eo)'./sum01',(sum00ee - sum00eo)'./sum00']); grid on;
% figure; plot(idle_times,[sum11ee'./sum11',sum10ee'./sum10',sum01ee'./sum01',sum00ee'./sum00']); grid on;
legend('R_{11}^{ee}','R_{10}^{ee}','R_{01}^{ee}','R_{00}^{ee}');


function [state, parity] = evolve_state(state, parity, t_idle, t_qp, correlation, T1, dt)
    for i = 1:t_idle/dt
        if rand() < dt/T1
            state = 0;
        elseif rand() < dt/t_qp
            parity = ~parity;
            if rand() < correlation
                state = ~state;
            end
        end
    end
end