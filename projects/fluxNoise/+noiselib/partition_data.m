function [seg_even, seg_odd] = partition_data(data)
%Partition data into even and odd for CPSD on single data set

seg_even = data(1:2:end);
seg_odd  = data(2:2:end);

end