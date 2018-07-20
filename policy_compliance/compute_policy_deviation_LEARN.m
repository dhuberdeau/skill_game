% Compute and save the policy deviation maps for beginning and asymptote
% windows.
%
% David Huberdeau 07/20/2018

init_wind = {1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050, 1201:1250, 1401:1450, 1601:1650, 1801:1850};
asym_wind = {151:200, 351:400, 551:600, 751:800, 951:1000, 1151:1200, 1351:1400, 1551:1600, 1751:1800, 1951:2000};

%% for successes:
[zu_trial_init_successes, zu_trial_asym_successes] = compute_policy_deviation_map(init_wind, asym_wind, 1);
save deviation_maps_LEARN_SUCCESSES zu_trial_init_successes zu_trial_asym_successes

%% for failures:
[zu_trial_init_failures, zu_trial_asym_failures] = compute_policy_deviation_map(init_wind, asym_wind, 0);
save deviation_maps_LEARN_FAILURES zu_trial_init_failures zu_trial_asym_failures