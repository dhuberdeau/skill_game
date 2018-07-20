% Compute and save the policy deviation maps for probe and pre-probe
% windows.
%
% David Huberdeau 07/20/2018

pre_wind = {26:75, 226:275, 426:475, 626:675, 826:875, 1026:1075, 1226:1275, 1426:1475, 1626:1675, 1826:1875};
probe_wind = {76:125, 276:325, 476:525, 676:725, 876:925, 1076:1125, 1276:1325, 1476:1525, 1676:1725, 1876:1925};

%% for successes:
[zu_trial_pre_successes, zu_trial_probe_successes] = compute_policy_deviation_map(pre_wind, probe_wind, 1);
save deviation_maps_PROBE_SUCCESSES zu_trial_pre_successes zu_trial_probe_successes

%% for failures:
[zu_trial_pre_failures, zu_trial_probe_failures] = compute_policy_deviation_map(pre_wind, probe_wind, 0);
save deviation_maps_PROBE_FAILURES zu_trial_pre_failures zu_trial_probe_failures
