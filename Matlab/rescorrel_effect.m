function [data] = rescorrel_effect(rc, a, b)
% Returns a plottable dataset of the "correlation" between each residue in a
% with all of those in b, according to rescorrel matrix rc

data = zeros(length(a),1);

for i = a
    data(i) = sum(rc(i,b));
end
