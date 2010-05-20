% Residue correlation factors for nucleic acids against the rest of the protein
rf_ttfm = rescorrel_factor(ttfm, 2:18, 34:718);
rf_tefm = rescorrel_factor(tefm, 2:13, 25:709);
rf_tec2a = rescorrel_factor(tec2a, 2:13, 25:709);
rf_teg3t = rescorrel_factor(teg3t, 2:13, 25:709);
rf_tea4c = rescorrel_factor(tea4c, 2:13, 25:709);
rf_tea4t = rescorrel_factor(tea4t, 2:13, 25:709);
rf_teg6c = rescorrel_factor(teg6c, 2:13, 25:709);
rf_tet7g = rescorrel_factor(tet7g, 2:13, 25:709);
rf_ttpg = rescorrel_factor(ttpg, 2:18, 19:703);
rf_aafm = rescorrel_factor(aafm, 704:724, 1:703);

%plot(1:17, rf_ttfm, 1:12, rf_tefm, 1:12, rf_tec2a, 1:12, rf_teg3t, 1:12, rf_tea4c, 1:12, rf_tea4t, 1:12, rf_teg6c, 1:12, rf_tet7g, 1:21, rf_aafm)
plot(1:17, rf_ttfm, 1:12, rf_tefm, 1:12, rf_tec2a, 1:12, rf_teg3t, 1:12, rf_tea4c, 1:12, rf_tea4t, 1:12, rf_teg6c, 1:12, rf_tet7g)
legend('3HK2 FM', '3F73 FM', 'C2A', 'G3T', 'A4C', 'A4T', 'G6C', 'T7G');
xlabel('Guide strand position');
ylabel('Correlation factor');

% Nucleic acid vs itself
nrf_ttfm = rescorrel_factor(ttfm, 2:18, 2:33);
nrf_tefm = rescorrel_factor(tefm, 2:13, 2:33);
nrf_tec2a = rescorrel_factor(tec2a, 2:13, 2:24);
nrf_teg3t = rescorrel_factor(teg3t, 2:13, 2:24);
nrf_tea4c = rescorrel_factor(tea4c, 2:13, 2:24);
nrf_tea4t = rescorrel_factor(tea4t, 2:13, 2:24);
nrf_teg6c = rescorrel_factor(teg6c, 2:13, 2:24);
nrf_tet7g = rescorrel_factor(tet7g, 2:13, 2:24);
nrf_ttpg = rescorrel_factor(ttpg, 2:18, 2:18);
nrf_aafm = rescorrel_factor(aafm, 704:724, 704:724);

%plot(1:17, nrf_ttfm, 1:12, nrf_tefm, 1:12, nrf_tec2a, 1:12, nrf_teg3t, 1:12, nrf_tea4c, 1:12, nrf_tea4t, 1:12, nrf_teg6c, 1:12, nrf_tet7g)
legend('3HK2 FM', '3F73 FM', 'C2A', 'G3T', 'A4C', 'A4T', 'G6C', 'T7G');
xlabel('Guide strand position');
ylabel('Correlation factor');
