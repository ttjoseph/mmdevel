function [data] = rescorrel_factor_distshell(rc, coords, a, b, inner, thickness)
% Returns a plottable dataset of the "correlation" between each residue in a
% with those in b that are between inner and inner+thickness away from the given residue in a
% according to rescorrel matrix rc and list of coordinates coords. Values are normalized by
% the number of residues taken into account

data = zeros(length(a),1);
dists = dist(coords');

for i = 1:length(a)
    b_in_shell = [];
    for j = 1:length(b)
        if dists(i, j) >= inner && dists(i, j) < (inner+thickness)
            b_in_shell = [b_in_shell j];
        end
    end
    data(i) = sum(rc(a(i), b_in_shell)) / length(b_in_shell);
end