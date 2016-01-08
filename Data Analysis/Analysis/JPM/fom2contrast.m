function contrast = fom2contrast(alpha)
% CONTRAST = fom2contrast(ALPHA) compute upper contrast boundry for given
% figure of merit value ALPHA.

contrast = alpha * (alpha + 1)^(-1 - 1 / alpha);

end