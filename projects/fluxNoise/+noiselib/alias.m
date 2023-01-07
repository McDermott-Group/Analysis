function x = alias(x0, bound)
% Takes a value or array x0 and aliases it into the range (-bound,bound]

x = mod(bound + x0, 2*bound) - bound;

end