% Load the data files - this is in a separate M file
% because loading those eigenvectors can take a long time

% Eigenvectors
ev_N = load('N.PCA.eigenvectors.dat');
ev_C = load('C.PCA.eigenvectors.dat');
ev_CA = load('CA.PCA.eigenvectors.dat');
orig_CB = load('CB.PCA.eigenvectors.dat');

% Coordinates of relevant atoms
solute_N = load('solute_N.pdb.dat');
solute_C = load('solute_C.pdb.dat');
solute_CA = load('solute_CA.pdb.dat');
% This one does not have NaNs for the coords of CBs that don't exist
orig_solute_CB = load('solute_CB.pdb.dat');

% Get rid of first column in solute_ which is the residue number.
% We can do this because we've previously ensured (outside Matlab)
% that the residue numbering has no gaps.
solute_N = solute_N(:, 2:4);
solute_C = solute_C(:, 2:4);
solute_CA = solute_CA(:, 2:4);
orig_solute_CB = orig_solute_CB(:, 2:4);

% Indices of missing CB atoms (due to glycines)
missing_CB = load('missing_CB.dat');
% Add zeros in empty spots due to CB atoms
fix_ev_CB;