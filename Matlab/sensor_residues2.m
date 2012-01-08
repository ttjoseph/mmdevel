% Alternative method for finding sensor residues in TtAgo

% Load files
fm = load('-ascii', 'MD/FullyMatched/rescorrel.txt'); fm = fm(25:709, 25:709);
g3t = load('-ascii', 'MD/G3T/rescorrel.txt'); g3t = g3t(25:709, 25:709);
a4c = load('-ascii', 'MD/A4C/rescorrel.txt'); a4c = a4c(25:709, 25:709);
a4t = load('-ascii', 'MD/A4T/rescorrel.txt'); a4t = a4t(25:709, 25:709);
g6c = load('-ascii', 'MD/G6C/rescorrel.txt'); g6c = g6c(25:709, 25:709);
t7g = load('-ascii', 'MD/T7G/rescorrel.txt'); t7g = t7g(25:709, 25:709);

% Calculate differences
NUM_MISMATCHES = 5;
dg3t = fm - g3t;
da4c = fm - a4c;
da4t = fm - a4t;
dg6c = fm - g6c;
dt7g = fm - t7g;

% Use sigma units
dg3t = (dg3t - mean(dg3t(:))) / std(dg3t(:));
da4c = (da4c - mean(da4c(:))) / std(da4c(:));
da4t = (da4t - mean(da4t(:))) / std(da4t(:));
dg6c = (dg6c - mean(dg6c(:))) / std(dg6c(:));
dt7g = (dt7g - mean(dt7g(:))) / std(dt7g(:));

% Descriptive statistics
dsum = (dg3t + da4c + da4t + dg6c + dt7g);
dsumsq = (dg3t.^2 + da4c.^2 + da4t.^2 + dg6c.^2 + dt7g.^2);
dmean = dsum / NUM_MISMATCHES;
% Variance for each cell in the matrix
dvar = (dsumsq - (dsum.*dmean)) / (NUM_MISMATCHES - 1);

% Root mean square (more informative than a straight sum) since opposite signs can't cancel
sumdrms = sum(sqrt(dsumsq / NUM_MISMATCHES));
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

