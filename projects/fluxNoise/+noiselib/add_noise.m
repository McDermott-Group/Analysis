function [es] = add_noise(es)

es = es + randn(length(es),1) * sqrt(3)*1e-2;

end