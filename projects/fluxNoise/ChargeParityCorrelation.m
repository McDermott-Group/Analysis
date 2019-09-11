data = loadMeasurementData;

% M1a = data.Single_Shot_Occupations_1(:,1:end-1);
% M1b = data.Single_Shot_Occupations_1(:,1:end-1);
% M2a = data.Single_Shot_Occupations_2(:,2:end);
% M2b = data.Single_Shot_Occupations_2(:,2:end);
% tau = data.Initialization_Time;
% n = length(data.Repetition_Index) - 1;

M1a = data.Single_Shot_Occupations_1;
M1b = data.Single_Shot_Occupations_2;
M2a = data.Single_Shot_Occupations_3;
M2b = data.Single_Shot_Occupations_4;
tau = data.Delay_Between_Readouts_2;
n = length(data.Repetition_Index);

P1 = mod(M1a + M1b, 2);
P2 = mod(M2a + M2b, 2);
dP = abs(P1-P2);

n11eo = zeros(length(tau),1);
n10ee = zeros(length(tau),1);
n10eo = zeros(length(tau),1);
n00eo = zeros(length(tau),1);
n01ee = zeros(length(tau),1);
n01eo = zeros(length(tau),1);

for i = 1:length(tau)
    for j = 1:n
        if (M1b(i,j) == 1 && M2a(i,j) == 1 && dP(i,j)==1)
            n11eo(i) = n11eo(i) + 1;
        elseif (M1b(i,j) == 1 && M2a(i,j) == 0 && dP(i,j)==0)
            n10ee(i) = n10ee(i) + 1;
        elseif (M1b(i,j) == 1 && M2a(i,j) == 0 && dP(i,j)==1)
            n10eo(i) = n10eo(i) + 1;
        elseif (M1b(i,j) == 0 && M2a(i,j) == 0 && dP(i,j)==1)
            n00eo(i) = n00eo(i) + 1;
        elseif (M1b(i,j) == 0 && M2a(i,j) == 1 && dP(i,j)==0)
            n01ee(i) = n01ee(i) + 1;
        elseif (M1b(i,j) == 0 && M2a(i,j) == 1 && dP(i,j)==1)
            n01eo(i) = n01eo(i) + 1;
        end
    end
end

figure();
plot(tau, n11eo, tau, n10ee, tau, n10eo, tau, n00eo, tau, n01ee, tau, n01eo)