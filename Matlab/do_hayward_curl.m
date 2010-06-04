function [the_curl] = do_hayward_curl(this_ev, ev_N, ev_C, ev_CA, ev_CB, solute_N, solute_C, solute_CA, solute_CB)
% Discrete curl - Tom Joseph <thomas.joseph@mssm.edu>
%
% This script implements the discrete curl method described in
% Model-free methods of analyzing domain motions in simulation...
% Hayward et al, Proteins 1997.
%
% Requires accompanying scripts load_data_files.m and fix_ev_CB.m
% (the latter takes care of missing CB atoms)
num_coords = size(ev_N, 1);
num_atoms = size(solute_CA, 1);
dVx = zeros(num_atoms, 3);
dVy = zeros(num_atoms, 3);
dVz = zeros(num_atoms, 3);
the_curl = zeros(num_atoms, 3);

% Column indices for these coordinates in solute_* matrices.
% There isn't actually a need for the X:Z indexing below
% but I haven't used Matlab in a while and this makes the intent
% of the statements clearer to me.
X = 1; Y = 2; Z = 3;

atom = 1; % Atom counter
coord = 1; % EV coordinate counter
while coord < num_coords   
    % All these are calculated with respect to CA
    % Each of these is the right hand side of Eq B3 in Hayward
    dx_CA = [ ev_C(coord, this_ev)
              ev_N(coord, this_ev)
              ev_CB(coord, this_ev) ] - ev_CA(coord, this_ev);
    
    dy_CA = [ ev_C(coord + 1, this_ev)
              ev_N(coord + 1, this_ev)
              ev_CB(coord + 1, this_ev) ] - ev_CA(coord + 1, this_ev);

    dz_CA = [ ev_C(coord + 2, this_ev)
              ev_N(coord + 2, this_ev)
              ev_CB(coord + 2, this_ev) ] - ev_CA(coord + 2, this_ev);
  
    % Now to make the matrix of differences in absolute coordinates
    % of the atoms in question
    % Order: C N CB in rows, x y z in cols
    xyz_CA = [ solute_CA(atom, X:Z)
               solute_CA(atom, X:Z)
               solute_CA(atom, X:Z) ];
    
    coord_diff_CA = [ solute_C(atom, X:Z)
                      solute_N(atom, X:Z)
                      solute_CB(atom, X:Z) ] - xyz_CA;
    
    % Record this set of partials if not any of the coord differences
    % are nonexistent (NaN, used for nonexistent CB atoms)
    % In English:
    % If this piece of the calculation doesn't include a nonexistent
    % CB atom, this set of partials is valid, so record it
    if ~any(isnan(coord_diff_CA))
        % dVx(atom, X:Z) = coord_diff_CA \ dx_CA;
        % dVy(atom, X:Z) = coord_diff_CA \ dy_CA;
        % dVz(atom, X:Z) = coord_diff_CA \ dz_CA;
        inv_coord_diff_CA = inv(coord_diff_CA);
        dVx(atom, X:Z) = inv_coord_diff_CA * dx_CA;
        dVy(atom, X:Z) = inv_coord_diff_CA * dy_CA;
        dVz(atom, X:Z) = inv_coord_diff_CA * dz_CA;
    end
    % Calculate the curl vector at this point
    the_curl(atom, X:Z) = [ dVz(atom, Y) - dVy(atom, Z)
                            dVx(atom, Z) - dVz(atom, X)
                            dVy(atom, X) - dVx(atom, Y) ];
    
    atom = atom + 1;
    coord = coord + 3;
end