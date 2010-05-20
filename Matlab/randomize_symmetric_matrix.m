function [m2] = randomize_symmetric_matrix(m)
%RANDOMIZE_SYMMETRIC_MATRIX   Swaps rows and columns of M at random, COUNT times.
%   This is used to generate a "random" matrix from an existing matrix
%   for bootstrapping purposes. Rows/columns are chosen at random from
%   the source matrix, with replacement.

if size(m, 1) ~= size(m, 2)
    error('Matrix must be square (and symmetric, but we don''t check for that).')
end

m2 = m;
dimension = size(m, 1);

for dest_idx = 1:dimension
    % Choose row/column index from source
    src_idx = floor(rand * dimension) + 1;
    % Copy over row
    m2(dest_idx, :) = m(src_idx, :);
    % Copy over column
    m2(:, dest_idx) = m(:, src_idx);
end

