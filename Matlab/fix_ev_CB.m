% Insert NaNs into ev_CB so that an index into ev_CB can be used for
% the other ev_* arrays
% Requires missing_CB and ev_CB to already be loaded from files
num_coords = size(orig_CB, 1) + size(missing_CB, 1) * 3;
num_orig_coords = size(orig_CB, 2);
ev_CB = zeros(num_coords, num_orig_coords);
num_atoms = num_coords / 3;
solute_CB = zeros(num_atoms, 3);
counter_CB = 1;
for atom = 1:num_atoms
    % We loop over the uniform atom indices
    % If it isn't one that isn't in ev_CB, copy
    dest_idx = (atom - 1) * 3 + 1 : (atom - 1 )* 3 + 3;
    src_idx = counter_CB : counter_CB + 2;
    if ~any(missing_CB == atom)
        %disp(sprintf('Saving atom %d : dest %d src %d', atom, dest_idx(1), src_idx(1)))
        ev_CB(dest_idx, :) = orig_CB(src_idx, :);
        solute_CB(atom, :) = orig_solute_CB((counter_CB - 1) / 3 + 1, :);
        counter_CB = counter_CB + 3;
    end
end
