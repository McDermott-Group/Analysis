function [es_jump, es_smooth] = filter_jumps(es,delta)

thresh = 0.07;%0.001;
es_jump = zeros(length(es),1) + es(1);
es_smooth = zeros(length(es),1) + es(1);
for i = 1:length(delta)
    if abs(delta(i)) < thresh
        es_smooth(i+1) = es_smooth(i) + delta(i);
        es_jump(i+1) = es_jump(i);
    else
        es_jump(i+1) = es_jump(i) + delta(i);
        es_smooth(i+1) = es_smooth(i);
    end
end

end