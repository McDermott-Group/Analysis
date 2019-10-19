function [seg_even, seg_odd] = chop_data(seg)
%Partition data into even and odd for CPSD on single data set

seg_even = seg(1:2:end);
seg_odd = seg(2:2:end);

end