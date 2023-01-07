function [ccf, lags] = crosscorrelate(x, y, n)
% Finds the cross correlation between the ones in two binary strings x and
% y of the same length.  n is the number of lags.  ccf is the cross
% correlated result, and lags is the lags=-n:n.  normalized to number of
% ones in x.

ccf = zeros(1,2*n+1);
ccf(n+1:end) = ccorr(x, y, n);
ccf(1:n+1) = flip(ccorr(y, x, n));
lags = -n:n;

end


function [acf] = ccorr(x, y, n)

acf = zeros(n+1,1);
for i = 0:n
    if sum(x(1:end-i)) == 0
        acf(i+1) = 0;
    else
        % we divide by the sum (not length) so that it is normalized to the
        % number of 1's and we get 100% corr if the 1s are correlated
        acf(i+1) = sum( x(1:end-i).*y(1+i:end) ) / sum(x(1:end-i));
    end
end

end