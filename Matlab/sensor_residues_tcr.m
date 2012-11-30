% Alternative method for finding sensor residues in TtAgo

% Load files
tcralone = load('-ascii', 'TCR_alone/rescorrel.txt');
taxwt = load('-ascii', 'TCR_TAXwt/rescorrel.txt');
p6a = load('-ascii', 'TCR_P6A/rescorrel.txt');
y5f = load('-ascii', 'TCR_Y5F/rescorrel.txt');
peptide_indices = 379:387;

% Calculate differences
NUM_STRUCTS = 3;
p6a_taxwt = p6a - taxwt;
y5f_taxwt = y5f - taxwt;
y5f_p6a = y5f - p6a;

% Use sigma units
p6a_taxwt = (p6a_taxwt - mean(p6a_taxwt(:))) / std(p6a_taxwt(:));
y5f_taxwt = (y5f_taxwt - mean(y5f_taxwt(:))) / std(y5f_taxwt(:));
y5f_p6a = (y5f_p6a - mean(y5f_p6a(:))) / std(y5f_p6a(:));

% Descriptive statistics
dsum = (p6a_taxwt + y5f_taxwt + y5f_p6a);
dsumsq = (p6a_taxwt.^2 + y5f_taxwt.^2 + y5f_p6a.^2);
dmean = dsum / NUM_STRUCTS;
% Variance for each cell in the matrix
dvar = (dsumsq - (dsum.*dmean)) / (NUM_STRUCTS - 1);

% Root mean square (more informative than a straight sum) since opposite signs can't cancel
sumdrms = sum(sqrt(dsumsq / NUM_STRUCTS));
sensor_residues = find(sumdrms > (mean(sumdrms) + 2*std(sumdrms)))

% How do we know which residues are relevant?


% NUM_SUBPLOTS = 5;
% BOUNDS = [1 685];
% RESIDUE = 546;
% figure;
% subplot(NUM_SUBPLOTS, 1, 1); 
% 
% plot(g3t(RESIDUE,:) - fm(RESIDUE,:));
% xlim(BOUNDS);
% ylabel('G3T');
% subplot(NUM_SUBPLOTS, 1, 2); 
% plot(a4c(RESIDUE,:) - fm(RESIDUE,:));
% xlim(BOUNDS);
% ylabel('A4C');
% subplot(NUM_SUBPLOTS, 1, 3); 
% plot(a4t(RESIDUE,:) - fm(RESIDUE,:));
% xlim(BOUNDS);
% ylabel('A4T');
% subplot(NUM_SUBPLOTS, 1, 4); 
% plot(g6c(RESIDUE,:) - fm(RESIDUE,:));
% xlim(BOUNDS);
% ylabel('G6C');
% subplot(NUM_SUBPLOTS, 1, 5); 
% plot(t7g(RESIDUE,:) - fm(RESIDUE,:));
% xlim(BOUNDS);
% ylabel('T7G');

