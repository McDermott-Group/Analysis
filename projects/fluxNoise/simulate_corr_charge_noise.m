pA = 0.05;
pB = 0.1;
pCorr = 0.03;

N = 1e7;

Corr = binornd(1, pCorr, 1, N);
A = Corr | binornd(1, pA, 1, N);
B = Corr | binornd(1, pB, 1, N);

nCorr = sum(A&B)
nCorrPredicted = N*(pCorr + 2*(pA-pCorr)*(pB-pCorr))