path(path,'~/devel/Matlab');
oldpwd = pwd(); % Save current directory
cd ~/TtAgoF/MD/FullyMatched
ttfm = load('-ascii', 'rescorrel.txt');
% Protein-only subset
ttfm_po = ttfm(34:718, 34:718);
cd ~/TtAgoF/MD/ProtGuide
ttpg = load('-ascii', 'rescorrel.txt');
ttpg_po = ttpg(19:703, 19:703);
cd ~/TtAgoF/MD/ProtOnly
ttpo = load('-ascii', 'rescorrel.txt');

% Truncated-nucleic-acid structures
cd ~/TtAgoE/MD/FullyMatched
tefm = load('-ascii', 'rescorrel.txt');
tefm_po = tefm(25:709, 25:709);
cd ~/TtAgoE/MD/ProtOnly
tepo = load('-ascii', 'rescorrel.txt');

cd ~/TtAgoE/MD/C2A
tec2a = load('-ascii', 'rescorrel.txt');
tec2a_po = tec2a(25:709, 25:709);
cd ~/TtAgoE/MD/G3T
teg3t = load('-ascii', 'rescorrel.txt');
teg3t_po = teg3t(25:709, 25:709);
cd ~/TtAgoE/MD/A4C
tea4c = load('-ascii', 'rescorrel.txt');
tea4c_po = tea4c(25:709, 25:709);
cd ~/TtAgoE/MD/A4T
tea4t = load('-ascii', 'rescorrel.txt');
tea4t_po = tea4t(25:709, 25:709);
cd ~/TtAgoE/MD/G6C
teg6c = load('-ascii', 'rescorrel.txt');
teg6c_po = teg6c(25:709, 25:709);
cd ~/TtAgoE/MD/T7G
tet7g = load('-ascii', 'rescorrel.txt');
tet7g_po = teg6c(25:709, 25:709);

% AaAgo structures
cd ~/AaAgo/ProtOnly
aapo = load('-ascii', 'rescorrel.txt');
cd ~/AaAgo/FullyMatchedRR
aafm = load('-ascii', 'rescorrel.txt');
aafm_po = aafm(1:703, 1:703);

% Get similarity value for comparing random matrices.
% We'll normalize to this.
%baseline = rescorrel_similarity(rand(685), rand(685));
baseline = 0.17; % Generated by bootstrapping

disp('Comparing protein-only subsets...')

% Calculate expected similarity values for randomly rearranged Ago rescorrel matrices.
% This is a bootstrapping-type method.
exv_ttfm_vs_ttpo = rescorrel_similarity_expected_value(ttfm_po, ttpo)
exv_ttfm_vs_tefm = rescorrel_similarity_expected_value(ttfm_po, tefm_po)
exv_ttfm_vs_ttpg = rescorrel_similarity_expected_value(ttfm_po, ttpg_po)

exv_ttpg_vs_ttpo = rescorrel_similarity_expected_value(ttpg_po, ttpo)
exv_tepo_vs_ttpo = rescorrel_similarity_expected_value(ttfm_po, tepo)
exv_tepo_vs_tefm = rescorrel_similarity_expected_value(tefm_po, tepo)

exv_tefm_vs_tec2a = rescorrel_similarity_expected_value(tefm_po, tec2a_po)
exv_tefm_vs_teg3t = rescorrel_similarity_expected_value(tefm_po, teg3t_po)
exv_tefm_vs_tea4c = rescorrel_similarity_expected_value(tefm_po, tea4c_po)
exv_tefm_vs_tea4t = rescorrel_similarity_expected_value(tefm_po, tea4t_po)
exv_tefm_vs_teg6c = rescorrel_similarity_expected_value(tefm_po, teg6c_po)
exv_tefm_vs_tet7g = rescorrel_similarity_expected_value(tefm_po, tet7g_po)

% Compare to TtAgoF protein-only simulation
ttfm_vs_ttpo = rescorrel_similarity(ttfm_po, ttpo)
ttfm_vs_tefm = rescorrel_similarity(ttfm_po, tefm_po)
ttfm_vs_ttpg = rescorrel_similarity(ttfm_po, ttpg_po)

ttpg_vs_ttpo = rescorrel_similarity(ttpg_po, ttpo)
tepo_vs_ttpo = rescorrel_similarity(ttfm_po, tepo)
tepo_vs_tefm = rescorrel_similarity(tefm_po, tepo)

tefm_vs_tec2a = rescorrel_similarity(tefm_po, tec2a_po)
tefm_vs_teg3t = rescorrel_similarity(tefm_po, teg3t_po)
tefm_vs_tea4c = rescorrel_similarity(tefm_po, tea4c_po)
tefm_vs_tea4t = rescorrel_similarity(tefm_po, tea4t_po)
tefm_vs_teg6c = rescorrel_similarity(tefm_po, teg6c_po)
tefm_vs_tet7g = rescorrel_similarity(tefm_po, tet7g_po)

% Compare AaAgo fully-matched to AaAgo protein only
aafm_vs_aapo = rescorrel_similarity(aafm_po, aapo)

cd(oldpwd); % Restore old working directory