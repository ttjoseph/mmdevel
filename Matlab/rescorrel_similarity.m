function [s] = rescorrel_similarity(rc1, rc2)
%RESCORREL_SIMILARITY   Calculate similarity or whatever.
%   If either of two corresponding cells have a value greater than the x-th percentile
%   value, then the other one must also exceed its x-th percentile value for the pair
%   to be considered similar. If neither cell has a value greater than threshold, the
%   pair is not considered.
%
%   Returns a fraction representing similarity, with a maximum of 1.
%   This result should be compared with the similaritiy of two random matrices,
%   by whatever measure of randomness you deem appropriate.
%
%   Tom Joseph <thomas.joseph@mssm.edu>

if size(rc1) ~= size(rc2)
    error('Sizes of rescorrel matrices must be the same.')
end

percentile = 0.926;

rc1v = sort(rc1(:));
rc2v = sort(rc2(:));
rc1_upper = rc1v(floor(length(rc1v) * percentile + 1));
rc2_upper = rc2v(floor(length(rc2v) * percentile + 1));

m = (rc1 > rc1_upper) & (rc2 > rc2_upper);
m2 = (rc1 > rc1_upper) | (rc2 > rc2_upper);
% imagesc(m2 - m)
s = sum(m(:)) / sum(m2(:));
