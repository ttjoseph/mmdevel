function [data] = rescorrel_factor(rc, a, b)
% Returns a plottable dataset of the "correlation" between each residue in a
% with all of those in b, according to rescorrel matrix rc

data = zeros(length(a),1);

for i = 1:length(a)
    data(i) = sum(rc(a(i),b));
end
