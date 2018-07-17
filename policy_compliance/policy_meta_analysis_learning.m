
score_fail_grp = cell(1,4);
score_succ_grp = cell(1,4);
slopes_pre_grp = cell(1, 4);
slopes_prb_grp = cell(1, 4);
score_fail_by_trial_grp = cell(2, 4); %row1: pre, row2: prb
score_succ_by_trial_grp = cell(2, 4);

for i_grp = 2:4
    
    pr_loc_gvn_succ_vLEARN_successes_v4;
    score_succ_grp{i_grp} = score;
%     score_fail_by_trial_grp{1, i_grp} = score_by_trial_pre;
%     score_fail_by_trial_grp{2, i_grp} = score_by_trial_prb;
    
    
    pr_loc_gvn_succ_vLEARN_failures_v4;
    score_fail_grp{i_grp} = score;
%     score_fail_by_trial_grp{1, i_grp} = score_by_trial_pre;
%     score_fail_by_trial_grp{2, i_grp} = score_by_trial_prb;
    
end

%%
score_fail_all = [[score_fail_grp{2}, nan(20, 10)];...
    [score_fail_grp{3}, nan(20, 10)];...
    score_fail_grp{4}];

score_succ_all = [[score_succ_grp{2}, nan(20, 10)];...
    [score_succ_grp{3}, nan(20, 10)];...
    score_succ_grp{4}];

%% plot
pre_wind_wind = [[.25 1.25, 2.25, 3.25 4.25], 5+[.25 1.25, 2.25, 3.25 4.25]];
prb_wind_wind = [[1 2, 3, 4 5], 5+[1 2, 3, 4 5]];
x = reshape([pre_wind_wind(1:grp_days(i_grp)); prb_wind_wind(1:grp_days(i_grp))], 2*grp_days(i_grp), 1);

figure;
subplot(2,1,1);
errorbar(x, nanmean(score_fail_all), sqrt(nanvar(score_fail_all)./sum(~isnan(score_fail_all))), 'k.-');
title('Failures')
axis([0 11 0.7 2.6])
subplot(2,1,2);
errorbar(x, nanmean(score_succ_all), sqrt(nanvar(score_succ_all)./sum(~isnan(score_succ_all))), 'k.-');
title('Successes')
axis([0 11 0.7 2.6])
%% export for analysis in R of linear mixed effect model

% within day analysis:
resp = [reshape(score_succ_all, 60*20, 1); reshape(score_fail_all, 60*20, 1)];
succ_fact = [ones(60*20, 1); zeros(60*20,1)];
wind_fact_continuous = repmat(reshape(repmat(1:20, 60, 1), 20*60, 1), 2, 1);
wind_fact_binary = repmat(reshape(repmat(repmat([1 2], 1, 10), 60, 1), 20*60, 1), 2, 1);
sub_fact = repmat(reshape(repmat((1:60)', 1, 20), 60*20, 1), 2, 1);
day_fact = repmat(reshape(repmat(reshape(repmat(1:10, 2, 1), 1, 20), 60, 1), 60*20, 1), 2, 1);

factors = [sub_fact, day_fact, wind_fact_binary, succ_fact]; 
response = resp;

csvwrite('policyDev_factors_LEARN_wDay', factors);
csvwrite('policyDev_response_LEARN_wDay', response);

% across day analysis:
resp = [reshape(score_succ_all(:, 2:19), 60*18, 1); reshape(score_fail_all(:, 2:19), 60*18, 1)];
succ_fact = [ones(60*18, 1); zeros(60*18,1)];
wind_fact_continuous = repmat(reshape(repmat(1:18, 60, 1), 18*60, 1), 2, 1);
wind_fact_binary = repmat(reshape(repmat(repmat([1 2], 1, 9), 60, 1), 18*60, 1), 2, 1);
sub_fact = repmat(reshape(repmat((1:60)', 1, 18), 60*18, 1), 2, 1);
day_fact = repmat(reshape(repmat(reshape(repmat(1:9, 2, 1), 1, 18), 60, 1), 60*18, 1), 2, 1);

factors = [sub_fact, day_fact, wind_fact_binary, succ_fact]; 
response = resp;

csvwrite('policyDev_factors_LEARN_aDay', factors);
csvwrite('policyDev_response_LEARN_aDay', response);
