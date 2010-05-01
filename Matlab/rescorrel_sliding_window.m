function [vals] = rescorrel_sliding_window(rescorrel)
% Sliding window average rescorrel value (really, a sliding box down the
% main diagonal...)

num_residues = size(rescorrel, 1);
halfwindowsize = 4;

vals = zeros(num_residues, 1);

for i = halfwindowsize + 1 : num_residues - halfwindowsize
    a = i - halfwindowsize;
    b = i + halfwindowsize - 1;
    vals(i) = mean(mean(rescorrel(a:b, a:b)));
end