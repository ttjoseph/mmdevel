function [o] = subspace_overlap(u, v)
% u and v and 3Natom x Nev dimensioned vectors
Nev = size(u, 2);
o = 0;
Nev = 50;

for i = 1:size(u,2)
    for j = 1:size(v,2)
        % uu = u(:, i) ./ norm(u(:, i));
        % vv = v(:, j) ./ norm(v(:, j));
        uu = u(:,i);
        vv = v(:,j);
        o = o + dot(uu, vv)^2;
    end
end