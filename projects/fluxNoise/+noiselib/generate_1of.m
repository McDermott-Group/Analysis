function [tdata] = generate_1of(reps, refreshTime, A, alpha)

dw = 1/(reps*refreshTime);
wdata = zeros(reps,1);
for i = 2:(length(wdata)/2)
    curr_val = A * (1/((i-1)*dw))^(alpha/2)*(normrnd(0,1) + 1i*normrnd(0,1));
%     phase = 2*pi*rand(1);
%     radius = normrnd(1,1);
%     curr_val = amplitude * (1/((i-1)*dw))^(alpha/2)*(radius*cos(phase) + radius*1i*sin(phase));
    wdata(i) = conj(curr_val);
    wdata(length(wdata) - i + 2) = curr_val;
end
wdata(length(wdata)/2+1) = A * (1/((i-1)*dw))^(alpha/2)*normrnd(0,1);
tdata = ifft(wdata);

end