function [delta] = calc_delta(es, stepsize)
    delta = es(stepsize+1:end) - es(1:end-stepsize);
end