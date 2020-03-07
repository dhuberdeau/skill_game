function [succ_meas_diff_set, kine_meas_diff_set] = test_policy_vs_succRate_interaction_b_fcn(score_succ_by_trial_grp, score_fail_by_trial_grp)
% During probe trials, distance travelled (and hazard) shows a well defined
% dip in performance. However, measures of kinematics (the variability and
% average trajectory, and the policy-deviation score) do not show the same
% reliable drop in metric. Why is this? One possibility is that kinematics
% haven't change (i.e. the ability to perform the learned movements is 
% retained during probes), but what has changed is the probability of a 
% failure. Another way to phrase this hypothesis is that we expect for 
% their to be an interaction between kinematic measures of performance and 
% measures of success. 
% However, testing for this interaction is not straighforward becasue the
% relationship between kinematics and task succcess are highly non-linear
% (because the track has an abrupt fall-off limit, etc.). So we need to
% some how establish the relationship, or, alternatively, subsample the
% kinematics among those that are nearly identical, and using that
% subsample again ask whether there has been a drop in performance. If
% there is truly an interaction between kinematics and success rate, the
% success rate for matched samples of kinematics should also show a
% difference between probe and pre-probe conditions. 
%
% This script does this analysis (for grps 3, 5, & 10, because grp 1 has 
% insufficient data for the analysis).
%
% David Huberdeau. 10/4/2017

load('exact_track_dist_full_v3.mat') %distance travelled measure
% load('score_by_trial_grp.mat') % policy deviation score

plot_tick_names = {'Kin', 'Succ'};

pre_probe_inds_by_grp = {26:75, 426:475, 826:875, 1826:1875};
probe_inds_by_grp = {76:125, 476:525, 876:925, 1876:1925};
probe_day_by_grp = {1, 3, 5, 10};

succ_meas = cell(2, 4);
for i_grp = 1:4
    succ_meas{1, i_grp} = exact_track_dist_fall{i_grp}(pre_probe_inds_by_grp{i_grp}, :);
    succ_meas{2, i_grp} = exact_track_dist_fall{i_grp}(probe_inds_by_grp{i_grp}, :);
end

kine_meas = cell(2, 4);
for i_grp = 1:4
    kine_meas{1, i_grp} = score_succ_by_trial_grp{1, i_grp}(:, :, probe_day_by_grp{i_grp});
    kine_meas{2, i_grp} = score_succ_by_trial_grp{2, i_grp}(:, :, probe_day_by_grp{i_grp});
end

n_subs = 0;
for i_grp = 1:4
    n_subs = n_subs + size(kine_meas{2,i_grp},2);
end

succ_meas_diff_set = cell(1,4);
kine_meas_diff_set = cell(1,4);
axis_min = 0;
axis_max = 0;

%%
succ_frac = nan(4, 21);
kin_frac = nan(4, 21);
for i_grp = 1:4
    succ_frac(i_grp, 1:size(succ_meas{1, i_grp},2)) = nanmean((succ_meas{2, i_grp})-nanmean(succ_meas{1, i_grp}))./nanmean(succ_meas{1, i_grp});
    kin_frac(i_grp, 1:size(succ_meas{1, i_grp},2)) = nanmean((kine_meas{2, i_grp})-nanmean(kine_meas{1, i_grp}))./nanmean(kine_meas{1, i_grp});
end
%% Match kinematics: what happens to distance travelled?
succ_meas_diff = nan(n_subs, 1);
kine_meas_diff = nan(n_subs, 1);
k_sub = 1;
TH_DIST = .05;
for i_grp = 2:4
    for i_sub = 1:size(kine_meas{1, i_grp},2)
        % collect Kinematic & success data for Probe window (called "Y")
        Y_temp = kine_meas{2, i_grp}(:, i_sub);
        Y = Y_temp(~isnan(Y_temp));
        succ_temp = succ_meas{2, i_grp}(:, i_sub);
        succ_Y = succ_temp(~isnan(Y_temp)); 
        if ~isempty(Y)
            % collect Kinematic & success data for Pre-probe window ("X")
            X_temp = kine_meas{1, i_grp}(:, i_sub);
            X = X_temp(~isnan(X_temp));
            [indx, distx] = knnsearch_unique(X, Y);

%             kine_meas_diff(k_sub) = nanmean(X(indx(distx > -TH_DIST & distx < TH_DIST))) - nanmean(Y(distx > -TH_DIST & distx < TH_DIST));
            kine_meas_diff(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
%             kine_meas_diff(k_sub) = nanmean(X(indx)) - nanmean(Y);
            
            succ_temp = succ_meas{1, i_grp}(:, i_sub);
            succ_X = succ_temp(~isnan(X_temp));
            succ_X_matched_samp = succ_X(indx(distx < TH_DIST));
            
            succ_meas_diff(k_sub) = nanmean(succ_X_matched_samp) - nanmean(succ_Y(distx < TH_DIST));
        end
        k_sub = 1 + k_sub;
    end
end

% disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% kine_meas_diff(disparate) = nan;
% succ_meas_diff(disparate) = nan;
[h_kin,p_kin,ci_kin,stat_kin] = ttest(kine_meas_diff)
[h_succ,p_succ,ci_succ,stat_succ] = ttest(succ_meas_diff)

succ_meas_diff_set{1} = succ_meas_diff;
kine_meas_diff_set{1} = kine_meas_diff;

axis_min = min([...
    nanmean(-[kine_meas_diff succ_meas_diff]) - sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_min]);
axis_max = max([...
    nanmean(-[kine_meas_diff succ_meas_diff]) + sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_max]);

%% Match distance travelled: what happens to kinematics
succ_meas_diff = nan(n_subs, 1);
kine_meas_diff = nan(n_subs, 1);
k_sub = 1;
TH_DIST = .01;
for i_grp = 2:4
    for i_sub = 1:size(kine_meas{1, i_grp},2)
        % collect Kinematic & success data for Probe window (called "Y")
        Y_temp = kine_meas{2, i_grp}(:, i_sub);
        Y = Y_temp(~isnan(Y_temp));
        succ_temp = succ_meas{2, i_grp}(:, i_sub);
        succ_Y = succ_temp(~isnan(Y_temp)); 
        if ~isempty(Y)
            % collect Kinematic & success data for Pre-probe window ("X")
            X_temp = kine_meas{1, i_grp}(:, i_sub);
            X = X_temp(~isnan(X_temp));
            
            succ_temp = succ_meas{1, i_grp}(:, i_sub);
            succ_X = succ_temp(~isnan(X_temp));
            
            [indx, distx] = knnsearch_unique(succ_X, succ_Y);

%             kine_meas_diff(k_sub) = nanmean(X(indx(distx > -TH_DIST & distx < TH_DIST))) - nanmean(Y(distx > -TH_DIST & distx < TH_DIST));
            kine_meas_diff(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
            succ_meas_diff(k_sub) = nanmean(succ_X(indx(distx < TH_DIST))) - nanmean(succ_Y(distx < TH_DIST));
        end
        k_sub = 1 + k_sub;
    end
end

% disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% kine_meas_diff(disparate) = nan;
% succ_meas_diff(disparate) = nan;
[h_kin,p_kin,ci_kin,stat_kin] = ttest(kine_meas_diff)
[h_succ,p_succ,ci_succ,stat_succ] = ttest(succ_meas_diff)

succ_meas_diff_set{2} = succ_meas_diff;
kine_meas_diff_set{2} = kine_meas_diff;

axis_min = min([...
    nanmean(-[kine_meas_diff succ_meas_diff]) - sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_min]);
axis_max = max([...
    nanmean(-[kine_meas_diff succ_meas_diff]) + sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_max]);
% 
% succ_meas_diff_2 = nan(n_subs, 1);
% kine_meas_diff_2 = nan(n_subs, 1);
% k_sub = 1;
% for i_grp = 2:4
%     for i_sub = 1:size(succ_meas{1, i_grp},2)
%         Y_temp = succ_meas{2, i_grp}(:, i_sub);
%         Y = Y_temp(~isnan(Y_temp));
%         kine_temp = kine_meas{2, i_grp}(:, i_sub);
%         kine_Y = kine_temp(~isnan(Y_temp));
%         if ~isempty(Y)
%             X_temp = succ_meas{1, i_grp}(:, i_sub);
%             X = X_temp(~isnan(X_temp));
% 
%             [indx, distx] = knnsearch_unique(X, Y);
% 
%             succ_meas_diff_2(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
% %             kine_meas_diff(k_sub) = nanmean(X(indx)) - nanmean(Y);
%             
%             kine_temp = kine_meas{1, i_grp}(:, i_sub);
%             kine_X = kine_temp(~isnan(X_temp));
%             kine_X_matched_samp = kine_X(indx(distx < TH_DIST));
%             
%             kine_meas_diff_2(k_sub) = nanmean(kine_X_matched_samp) - nanmean(kine_Y(distx < TH_DIST));
%         end
%         k_sub = 1 + k_sub;
%     end
% end
% 
% % disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% % kine_meas_diff(disparate) = nan;
% % succ_meas_diff(disparate) = nan;
% [h_kin_2,p_kin_2,ci_kin_2,stat_kin_2] = ttest(kine_meas_diff_2)
% [h_succ_2,p_succ_2,ci_succ_2,stat_succ_2] = ttest(succ_meas_diff_2)
% 
% 
% succ_meas_diff_set{2} = succ_meas_diff_2;
% kine_meas_diff_set{2} = kine_meas_diff_2;
% 
% axis_min = min([...
%     nanmean(-[kine_meas_diff_2 succ_meas_diff_2]) - sqrt(nanvar([kine_meas_diff_2 succ_meas_diff_2])./size(kine_meas_diff_2,1)),...
%     axis_min]);
% axis_max = max([...
%     nanmean(-[kine_meas_diff_2 succ_meas_diff_2]) + sqrt(nanvar([kine_meas_diff_2 succ_meas_diff_2])./size(kine_meas_diff_2,1)),...
%     axis_max]);

%%
 %%%%% FOR FAILURES:
 
 succ_meas = cell(2, 4);
for i_grp = 1:4
    succ_meas{1, i_grp} = exact_track_dist_fall{i_grp}(pre_probe_inds_by_grp{i_grp}, :);
    succ_meas{2, i_grp} = exact_track_dist_fall{i_grp}(probe_inds_by_grp{i_grp}, :);
end

kine_meas = cell(2, 4);
for i_grp = 1:4
    kine_meas{1, i_grp} = score_fail_by_trial_grp{1, i_grp}(:, :, probe_day_by_grp{i_grp});
    kine_meas{2, i_grp} = score_fail_by_trial_grp{2, i_grp}(:, :, probe_day_by_grp{i_grp});
end
 
 %% Match kinematics: what happens to distance travelled?
succ_meas_diff = nan(n_subs, 1);
kine_meas_diff = nan(n_subs, 1);
k_sub = 1;
for i_grp = 2:4
    for i_sub = 1:size(kine_meas{1, i_grp},2)
        % collect Kinematic & success data for Probe window (called "Y")
        Y_temp = kine_meas{2, i_grp}(:, i_sub);
        Y = Y_temp(~isnan(Y_temp));
        succ_temp = succ_meas{2, i_grp}(:, i_sub);
        succ_Y = succ_temp(~isnan(Y_temp)); 
        if ~isempty(Y)
            % collect Kinematic & success data for Pre-probe window ("X")
            X_temp = kine_meas{1, i_grp}(:, i_sub);
            X = X_temp(~isnan(X_temp));
            [indx, distx] = knnsearch_unique(X, Y);

%             kine_meas_diff(k_sub) = nanmean(X(indx(distx > -TH_DIST & distx < TH_DIST))) - nanmean(Y(distx > -TH_DIST & distx < TH_DIST));
            kine_meas_diff(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
%             kine_meas_diff(k_sub) = nanmean(X(indx)) - nanmean(Y);
            
            succ_temp = succ_meas{1, i_grp}(:, i_sub);
            succ_X = succ_temp(~isnan(X_temp));
            succ_X_matched_samp = succ_X(indx(distx < TH_DIST));
            
            succ_meas_diff(k_sub) = nanmean(succ_X_matched_samp) - nanmean(succ_Y(distx < TH_DIST));
        end
        k_sub = 1 + k_sub;
    end
end

% disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% kine_meas_diff(disparate) = nan;
% succ_meas_diff(disparate) = nan;
[h_kin,p_kin,ci_kin,stat_kin] = ttest(kine_meas_diff)
[h_succ,p_succ,ci_succ,stat_succ] = ttest(succ_meas_diff)

succ_meas_diff_set{3} = succ_meas_diff;
kine_meas_diff_set{3} = kine_meas_diff;

axis_min = min([...
    nanmean(-[kine_meas_diff succ_meas_diff]) - sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_min]);
axis_max = max([...
    nanmean(-[kine_meas_diff succ_meas_diff]) + sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_max]);

%% Match distance travelled: what happens to kinematics
succ_meas_diff = nan(n_subs, 1);
kine_meas_diff = nan(n_subs, 1);
k_sub = 1;
TH_DIST = .01;
for i_grp = 2:4
    for i_sub = 1:size(kine_meas{1, i_grp},2)
        % collect Kinematic & success data for Probe window (called "Y")
        Y_temp = kine_meas{2, i_grp}(:, i_sub);
        Y = Y_temp(~isnan(Y_temp));
        succ_temp = succ_meas{2, i_grp}(:, i_sub);
        succ_Y = succ_temp(~isnan(Y_temp)); 
        if ~isempty(Y)
            % collect Kinematic & success data for Pre-probe window ("X")
            X_temp = kine_meas{1, i_grp}(:, i_sub);
            X = X_temp(~isnan(X_temp));
            
            succ_temp = succ_meas{1, i_grp}(:, i_sub);
            succ_X = succ_temp(~isnan(X_temp));
            
            [indx, distx] = knnsearch_unique(succ_X, succ_Y);

%             kine_meas_diff(k_sub) = nanmean(X(indx(distx > -TH_DIST & distx < TH_DIST))) - nanmean(Y(distx > -TH_DIST & distx < TH_DIST));
            kine_meas_diff(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
            succ_meas_diff(k_sub) = nanmean(succ_X(indx(distx < TH_DIST))) - nanmean(succ_Y(distx < TH_DIST));
        end
        k_sub = 1 + k_sub;
    end
end

% disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% kine_meas_diff(disparate) = nan;
% succ_meas_diff(disparate) = nan;
[h_kin,p_kin,ci_kin,stat_kin] = ttest(kine_meas_diff)
[h_succ,p_succ,ci_succ,stat_succ] = ttest(succ_meas_diff)

succ_meas_diff_set{4} = succ_meas_diff;
kine_meas_diff_set{4} = kine_meas_diff;

axis_min = min([...
    nanmean(-[kine_meas_diff succ_meas_diff]) - sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_min]);
axis_max = max([...
    nanmean(-[kine_meas_diff succ_meas_diff]) + sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
    axis_max]);
%%
% 
% succ_meas = cell(2, 4);
% for i_grp = 1:4
%     succ_meas{1, i_grp} = exact_track_dist_fall{i_grp}(pre_probe_inds_by_grp{i_grp}, :);
%     succ_meas{2, i_grp} = exact_track_dist_fall{i_grp}(probe_inds_by_grp{i_grp}, :);
% end
% 
% kine_meas = cell(2, 4);
% for i_grp = 1:4
%     kine_meas{1, i_grp} = score_fail_by_trial_grp{1, i_grp}(:, :, probe_day_by_grp{i_grp});
%     kine_meas{2, i_grp} = score_fail_by_trial_grp{2, i_grp}(:, :, probe_day_by_grp{i_grp});
% end
% 
% n_subs = 0;
% for i_grp = 1:4
%     n_subs = n_subs + size(kine_meas{2,i_grp},2);
% end
% 
% %%
% succ_frac = nan(4, 21);
% kin_frac = nan(4, 21);
% for i_grp = 1:4
%     succ_frac(i_grp, 1:size(succ_meas{1, i_grp},2)) = nanmean((succ_meas{2, i_grp})-nanmean(succ_meas{1, i_grp}))./nanmean(succ_meas{1, i_grp});
%     kin_frac(i_grp, 1:size(succ_meas{1, i_grp},2)) = nanmean((kine_meas{2, i_grp})-nanmean(kine_meas{1, i_grp}))./nanmean(kine_meas{1, i_grp});
% end
% %% Match kinematics: what happens to distance travelled?
% succ_meas_diff = nan(n_subs, 1);
% kine_meas_diff = nan(n_subs, 1);
% k_sub = 1;
% for i_grp = 2:4
%     for i_sub = 1:size(kine_meas{1, i_grp},2)
%         Y_temp = kine_meas{2, i_grp}(:, i_sub);
%         Y = Y_temp(~isnan(Y_temp));
%         succ_temp = succ_meas{2, i_grp}(:, i_sub);
%         succ_Y = succ_temp(~isnan(Y_temp));
%         if ~isempty(Y)
%             X_temp = kine_meas{1, i_grp}(:, i_sub);
%             X = X_temp(~isnan(X_temp));
% %             [indx, distx] = knnsearch(X, Y);
%             [indx, distx] = knnsearch_unique(X, Y);
%             kine_meas_diff(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
% %             kine_meas_diff(k_sub) = nanmean(X(indx)) - nanmean(Y);
%             
%             succ_temp = succ_meas{1, i_grp}(:, i_sub);
%             succ_X = succ_temp(~isnan(X_temp));
%             succ_X_matched_samp = succ_X(indx(distx < TH_DIST));
%             
%             succ_meas_diff(k_sub) = nanmean(succ_X_matched_samp) - nanmean(succ_Y(distx < TH_DIST));
%         end
%         k_sub = 1 + k_sub;
%     end
% end
% 
% % disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% % kine_meas_diff(disparate) = nan;
% % succ_meas_diff(disparate) = nan;
% [h_kin,p_kin,ci_kin,stat_kin] = ttest(kine_meas_diff)
% [h_succ,p_succ,ci_succ,stat_succ] = ttest(succ_meas_diff)
% 
% 
% 
% succ_meas_diff_set{3} = succ_meas_diff;
% kine_meas_diff_set{3} = kine_meas_diff;
% 
% axis_min = min([...
%     nanmean(-[kine_meas_diff succ_meas_diff]) - sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
%     axis_min]);
% axis_max = max([...
%     nanmean(-[kine_meas_diff succ_meas_diff]) + sqrt(nanvar([kine_meas_diff succ_meas_diff])./size(kine_meas_diff,1)),...
%     axis_max]);
% 
% %% Match distance travelled: what happens to kinematics
% succ_meas_diff_2 = nan(n_subs, 1);
% kine_meas_diff_2 = nan(n_subs, 1);
% k_sub = 1;
% for i_grp = 2:4
%     for i_sub = 1:size(succ_meas{1, i_grp},2)
%         Y_temp = succ_meas{2, i_grp}(:, i_sub);
%         Y = Y_temp(~isnan(Y_temp));
%         kine_temp = kine_meas{2, i_grp}(:, i_sub);
%         kine_Y = kine_temp(~isnan(Y_temp));
%         if ~isempty(Y)
%             X_temp = succ_meas{1, i_grp}(:, i_sub);
%             X = X_temp(~isnan(X_temp));
% %             [indx, distx] = knnsearch(X, Y);
%             [indx, distx] = knnsearch_unique(X, Y);
%             succ_meas_diff_2(k_sub) = nanmean(X(indx(distx < TH_DIST))) - nanmean(Y(distx < TH_DIST));
% %             kine_meas_diff(k_sub) = nanmean(X(indx)) - nanmean(Y);
%             
%             kine_temp = kine_meas{1, i_grp}(:, i_sub);
%             kine_X = kine_temp(~isnan(X_temp));
%             kine_X_matched_samp = kine_X(indx(distx < TH_DIST));
%             
%             kine_meas_diff_2(k_sub) = nanmean(kine_X_matched_samp) - nanmean(kine_Y(distx < TH_DIST));
%         end
%         k_sub = 1 + k_sub;
%     end
% end
% 
% % disparate = kine_meas_diff > .2 | kine_meas_diff < -.2;
% % kine_meas_diff(disparate) = nan;
% % succ_meas_diff(disparate) = nan;
% [h_kin_2,p_kin_2,ci_kin_2,stat_kin_2] = ttest(kine_meas_diff_2)
% [h_succ_2,p_succ_2,ci_succ_2,stat_succ_2] = ttest(succ_meas_diff_2)
% 
% 
% succ_meas_diff_set{4} = succ_meas_diff_2;
% kine_meas_diff_set{4} = kine_meas_diff_2;
% 
% axis_min = min([...
%     nanmean(-[kine_meas_diff_2 succ_meas_diff_2]) - sqrt(nanvar([kine_meas_diff_2 succ_meas_diff_2])./size(kine_meas_diff_2,1)),...
%     axis_min]);
% axis_max = max([...
%     nanmean(-[kine_meas_diff_2 succ_meas_diff_2]) + sqrt(nanvar([kine_meas_diff_2 succ_meas_diff_2])./size(kine_meas_diff_2,1)),...
%     axis_max]);

%% do all plotting:

figure; ah = axes; 
for i_type = 1:4
    subplot(2,2,i_type);
    hold on;
    kine_meas_diff = kine_meas_diff_set{i_type};
    succ_meas_diff = succ_meas_diff_set{i_type};
    errorbar([1,2], nanmean(-[kine_meas_diff succ_meas_diff]), sqrt(nanvar([kine_meas_diff succ_meas_diff])./[size(kine_meas_diff,1), size(succ_meas_diff,1)]),'.')
    bar([1,2], nanmean(-[kine_meas_diff succ_meas_diff]),'k')
    switch i_type
        case 1
            title('Matched Kinematics');
        case 2
            title('Matched Success');
        case 3
            title('Matched Kinematics');
        case 4
            title('Matched Success');
        otherwise
            warning('Invalid type specified for plotting');
    end
    axis([.5 2.5 axis_min axis_max])
    set(gca, 'xtick', 1:2, 'xticklabel', plot_tick_names);
end