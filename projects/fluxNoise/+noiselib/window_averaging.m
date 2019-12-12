function [apsd] = window_averaging(psd)
% apsd = psd;
%Averaging with f window

apsd = zeros(size(psd));
for i = 1:size(psd,1)
    filter_fl = max(1,round(i - i/4));
    filter_fh = min(size(psd,1),round(i + i/4));
    apsd(i) = mean(psd(filter_fl:filter_fh));
end

end