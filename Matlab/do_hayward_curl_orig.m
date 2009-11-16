% This script implements the discrete curl method described in
% Model-free methods of analyzing domain motions in simulation...
% Hayward et al, Proteins 1997.
num_coords = size(ev, 1);
this_ev = 1;
num_atoms = size(solute_CA, 1);
dVx = zeros(num_atoms, 3);
dVy = zeros(num_atoms, 3);
dVz = zeros(num_atoms, 3);
the_curl = zeros(num_atoms, 3);

% column indices for these coordinates in solute_* matrices
X = 1; Y = 2; Z = 3;

atom = 1; % Atom counter
i = 1; % EV coordinate counter
while i < num_coords
    
    while find(missing_CB == atom)
        i = i + 9;
        atom = atom + 1;
    end
    
    Nx = i; Ny = i + 1; Nz = i + 2;
    Cx = i + 3; Cy = i + 4; Cz = i + 5;
    CAx = i + 6; CAy = i + 7; CAz = i + 8;
    CBx = i + 9; CBy = i + 10; CBz = i + 11;
    % ev(N, 1) refers to the displacement of the ith N atom 
    % in the first eigenvector
    
    % All these are calculated with respect to N
    % Each of these is the right hand side of Eq B3 in Hayward
    dx_N = [ ev(Cx, this_ev)
             ev(CAx, this_ev)
             ev(CBx, this_ev) ] - ev(Nx, this_ev);
    
    dy_N = [ ev(Cy, this_ev)
             ev(CAy, this_ev)
             ev(CBy, this_ev) ] - ev(Ny, this_ev);

    dz_N = [ ev(Cz, this_ev)
             ev(CAz, this_ev)
             ev(CBz, this_ev) ] - ev(Nz, this_ev);
  
    % Now to make the matrix of differences in absolute coordinates
    % of the atoms in question
    % Order: C CA CB in rows, x y z in cols
    xyz_N = [ solute_N(atom, X:Z)
              solute_N(atom, X:Z)
              solute_N(atom, X:Z) ];
    
    coord_diff_N = [ solute_C(atom, X:Z)
                     solute_CA(atom, X:Z)
                     solute_CB(atom, X:Z) ] - xyz_N;
  
    % Record this set of partials
    dVx(atom, X:Z) = dx_N \ coord_diff_N;
    dVy(atom, X:Z) = dy_N \ coord_diff_N;
    dVz(atom, X:Z) = dz_N \ coord_diff_N;
    
    % Calculate the curl vector at this point
    the_curl(atom, X:Z) = [ dVz(atom, Y) - dVy(atom, Z)
                            dVx(atom, Z) - dVz(atom, X)
                            dVy(atom, X) - dVx(atom, Y) ];
    
    atom = atom + 1;
    i = i + (4 * 3);
end