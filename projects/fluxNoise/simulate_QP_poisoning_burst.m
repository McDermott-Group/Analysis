T = 100e-6;
T1 = 25e-6;
gamma_up = 0.03;
gamma_down = exp(-T/T1);

N = 4000*10;

A = 20*gamma_up;
tau = 10e-3;
gamma = zeros(1,N);
gamma(5000:6000) = A*exp(-T*(0:1000)/tau);

r = rand(1,N);
o = zeros(1,N);
for i=2:N
    if o(i-1) == 0
        o(i) = (r(i) < gamma_up + gamma(i));
    elseif o(i-1) == 1
        o(i) = (r(i) < gamma_down + gamma(i));
    end
end

o = reshape(o,[4000,10]);

figure; imagesc(o)
figure; imagesc(movmean(o,10))

figure; plot(sum(o)/4000)

figure; hold on;
ccf_baseline = 1/9*noiselib.crosscorrelate(o(:,1),o(:,1),500);
for i = 3:10
    ccf_baseline = ccf_baseline + 1/9*noiselib.crosscorrelate(o(:,i),o(:,i),500);
end
ccf = noiselib.crosscorrelate(o(:,2),o(:,2),500);
plot(-500:500,ccf_baseline);
plot(-500:500, ccf);