% main_policy_analysis.m
%
% Script for policy analysis. Computes policy deviation for each group
% across available days.
%
% David Huberdeau, 12/1/2018

%% Load in data:
load('tilt_dir7.mat');
load('tilt_mag7.mat');
load('traj_dir7.mat');
load('traj_loc7.mat');
load('exact_track_dist_full_v3.mat');

%% define parameters:

% analysis window and subset of kinematic data:
len_range = [.25 .7]; % region of analysis more broadly defined.
% len_range = [.25 .45]; % do half the track length
% len_range = [.306 .4]; % region of analysis narrowly defined.

% note: tilt_dir, tilt_mag, traj_dir, and traj_loc were computed assuming
% the longest considered len_range window, .25 to .7. If you want to
% analyze a shorter window (e.g. to focus on the portion of the track with
% the highest hazard rate) will need to select sub-set of trajectory data:

range_inds = 1:size(traj_loc{1},1);
min_range = .25; max_range = .7;
range_resolution = (max_range - min_range)/length(range_inds);
range_vals = (min_range+range_resolution):range_resolution:max_range;

score_signal_window = 1:length(range_vals);
% score_signal_window = (length(range_vals) - 10):length(range_vals);

[~, sub_min_range] = min(abs(range_vals - len_range(1)));
[~, sub_max_range] = min(abs(range_vals - len_range(2)));
sub_range_inds = range_inds((sub_min_range):sub_max_range);

for i_grp = 1:length(tilt_dir)
    tilt_dir{i_grp} = tilt_dir{i_grp}(sub_range_inds, :, :);
    tilt_mag{i_grp} = tilt_mag{i_grp}(sub_range_inds, :, :);
    traj_dir{i_grp} = traj_dir{i_grp}(sub_range_inds, :, :);
    traj_loc{i_grp} = traj_loc{i_grp}(sub_range_inds, :, :);
end


ma_wind_size = 25;
group_colors = {'y', 'r', 'b', 'g'};
% probe_winds_indices = {4:5, 20:21, 36:37, 76:77};
probe_winds_indices = {5, 21, 37, 77};
axis_ranges = [0 2000 .5 1.4];
self_scale = 0;

% outpath = ['/Users/david/OneDrive/Documents/JHU/BLAM_lab/',...
%     'Projects/skillLearning_antGame_A/Matlab/public_analysis/data/'];
outpath = pwd;

% response_transform_func = @(x) log(log(10*x));
response_transform_func = @(x) x;
%% compute policy and scores:
policy = create_policy_map_L1O_fcn_v2(tilt_dir, tilt_mag, traj_dir, traj_loc, len_range);

%%
score_s = compute_policy_deviation_map_fcn(policy, exact_track_dist_fall, len_range, 1, ...
    tilt_dir, traj_dir, tilt_mag, traj_loc);
score_f = compute_policy_deviation_map_fcn(policy, exact_track_dist_fall, len_range, 0, ...
    tilt_dir, traj_dir, tilt_mag, traj_loc);
score_c = compute_policy_deviation_map_fcn(policy, exact_track_dist_fall, len_range, 2, ...
    tilt_dir, traj_dir, tilt_mag, traj_loc);

%% plot results:

fig_pos = [1, 796, 1000, 159];
bars_pos = [1 507 196 215];

training_winds_grp = {[1 1 1 0 0 1 1 1],...
    [ones(1, 19), 0, 0, ones(1, 19)],...
    [ones(1, 35), 0, 0, ones(1, 3)],...
    [ones(1, 40 + 35), 0, 0, ones(1, 3)]};

% successes:
successes_probe_diff_grp = nan(21, 4);
successes_avg_policy = nan(6000, 4);
k_ind = 1;

figure;
for i_grp = 1:4
    subplot(1,4,i_grp);
    hold on;
    
    grp_score = reshape(...
        nanmean(score_s{i_grp}(:,score_signal_window,:), 2),...
        size(score_s{i_grp},1),...
        size(score_s{i_grp},3));
    
    ma_score_ = mafilt(grp_score, ma_wind_size,1);
    ma_score = response_transform_func(exclude_outliers(ma_score_));
    
    % collect window data
    temp_winds = repmat(ma_wind_size*(1:size(ma_score,1)), 1, size(ma_score,2));
    temp_subj = repmat((1:size(ma_score,2)), size(ma_score, 1), 1);
    successes_avg_policy(k_ind:(k_ind - 1 + numel(ma_score)), :) = ...
        [ma_score(:), temp_winds(:),...
        i_grp*ones(numel(ma_score), 1),...
        temp_subj(:)];
    k_ind = k_ind + numel(ma_score);
    
    % collect probe diff data
    inds = 1:size(ma_score,1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;
    probe_winds = nanmean(ma_score(probe_inds, :));
    pre_p_winds = nanmean(ma_score(pre_p_inds, :));
    successes_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
    
    errorfield(...
        ma_wind_size*(1:size(ma_score,1)), ...
        nanmean(ma_score, 2), ...
        sqrt(nanvar(ma_score, [], 2)./sum(~isnan(ma_score),2)), ...
        group_colors{i_grp});
    
    errorfield(...
        ma_wind_size*probe_winds_indices{i_grp}, ...
        nanmean(ma_score(probe_winds_indices{i_grp},:), 2), ...
        sqrt(nanvar(ma_score(probe_winds_indices{i_grp},:), [], 2)./...
        sum(~isnan(ma_score(probe_winds_indices{i_grp},:)),2)), ...
        'k');
    title(['Group: ', num2str(i_grp), ', Succ']);
    if ~self_scale
        axis(axis_ranges)
    end
end
f2 = gcf;
set(f2, 'Position', fig_pos);
successes_avg_policy =...
    successes_avg_policy(~isnan(successes_avg_policy(:, 1)), :);
csvwrite([outpath, 'policy_successes'], successes_avg_policy);

diff_response_y = successes_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(successes_probe_diff_grp,1), 1), numel(successes_probe_diff_grp), 1), ...
    nan(numel(successes_probe_diff_grp), 1), ...
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])'];
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'policy_successes_probe_diff_grp'], [diff_response_y, diff_design_X]);

figure; hold on;
errorbar([1 2 3 4], nanmean(successes_probe_diff_grp), ...
    sqrt(nanvar(successes_probe_diff_grp)./sum(~isnan(successes_probe_diff_grp))), 'k.')
bar([1 2 3 4], nanmean(successes_probe_diff_grp));
title('successes')
axis([0 5 -.15 .1])
f2 = gcf;
set(f2, 'Position', bars_pos);

% failures:
failures_avg_policy = nan(6000, 4);
failures_probe_diff_grp = nan(21, 4);
axis_ranges = [0 2000 1 1.9];
k_ind = 1;
figure;
for i_grp = 1:4
    subplot(1,4,i_grp);
    hold on;
    
    grp_score = reshape(...
        nanmean(score_f{i_grp}(:,score_signal_window,:), 2),...
        size(score_f{i_grp},1),...
        size(score_f{i_grp},3));
    
    ma_score_ = mafilt(grp_score, ma_wind_size, 1);
    ma_score = response_transform_func(exclude_outliers(ma_score_));
    
    temp_winds = repmat(ma_wind_size*(1:size(ma_score,1)), 1, size(ma_score,2));
    temp_subj = repmat((1:size(ma_score,2)), size(ma_score, 1), 1);
    failures_avg_policy(k_ind:(k_ind - 1 + numel(ma_score)), :) = ...
        [ma_score(:), temp_winds(:),...
        i_grp*ones(numel(ma_score), 1),...
        temp_subj(:)];
    k_ind = k_ind + numel(ma_score);
    
    % collect probe diff data
    inds = 1:size(ma_score,1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;
    probe_winds = nanmean(ma_score(probe_inds, :));
    pre_p_winds = nanmean(ma_score(pre_p_inds, :));
    failures_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
    
    errorfield(...
        ma_wind_size*(1:size(ma_score,1)), ...
        nanmean(ma_score, 2), ...
        sqrt(nanvar(ma_score, [], 2)./sum(~isnan(ma_score),2)), ...
        group_colors{i_grp});
    
    errorfield(...
        ma_wind_size*probe_winds_indices{i_grp}, ...
        nanmean(ma_score(probe_winds_indices{i_grp},:), 2), ...
        sqrt(nanvar(ma_score(probe_winds_indices{i_grp},:), [], 2)./...
        sum(~isnan(ma_score(probe_winds_indices{i_grp},:)),2)), ...
        'k');
    title(['Group: ', num2str(i_grp), ', Fail']);
    if ~self_scale
        axis(axis_ranges)
    end
end
f1 = gcf;
set(f1, 'Position', fig_pos);
failures_avg_policy =...
    failures_avg_policy(~isnan(failures_avg_policy(:, 1)), :);
csvwrite([outpath, 'policy_failures'], failures_avg_policy);

diff_response_y = failures_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(failures_probe_diff_grp,1), 1), numel(failures_probe_diff_grp), 1), ...
    nan(numel(failures_probe_diff_grp), 1), ...
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])'];
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'policy_failures_probe_diff_grp'], [diff_response_y, diff_design_X]);

figure; hold on;
errorbar([1 2 3 4], nanmean(failures_probe_diff_grp), ...
    sqrt(nanvar(failures_probe_diff_grp)./sum(~isnan(failures_probe_diff_grp))), 'k.')
bar([1 2 3 4], nanmean(failures_probe_diff_grp));
title('failures')
axis([0 5 -.15 .1])
f1 = gcf;
set(f1, 'Position', bars_pos);

% combined:
union_avg_policy = nan(6000, 4);
union_probe_diff_grp = nan(21, 4);
k_ind = 1;
figure;
for i_grp = 1:4
    subplot(1,4,i_grp);
    hold on;
    
    grp_score = reshape(...
        nanmean(score_c{i_grp}(:,score_signal_window,:), 2),...
        size(score_c{i_grp},1),...
        size(score_c{i_grp},3));
    
    ma_score_ = mafilt(grp_score, ma_wind_size, 1);
    ma_score = response_transform_func(exclude_outliers(ma_score_));
    
    temp_winds = repmat(ma_wind_size*(1:size(ma_score,1)), 1, size(ma_score,2));
    temp_subj = repmat((1:size(ma_score,2)), size(ma_score, 1), 1);
    union_avg_policy(k_ind:(k_ind - 1 + numel(ma_score)), :) = ...
        [ma_score(:), temp_winds(:),...
        i_grp*ones(numel(ma_score), 1),...
        temp_subj(:)];
    k_ind = k_ind + numel(ma_score);
    
    % collect probe diff data
    inds = 1:size(ma_score,1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;
    probe_winds = nanmean(ma_score(probe_inds, :));
    pre_p_winds = nanmean(ma_score(pre_p_inds, :));
    union_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
    
    errorfield(...
        ma_wind_size*(1:size(ma_score,1)), ...
        nanmean(ma_score, 2), ...
        sqrt(nanvar(ma_score, [], 2)./sum(~isnan(ma_score),2)), ...
        group_colors{i_grp});
    
    errorfield(...
        ma_wind_size*probe_winds_indices{i_grp}, ...
        nanmean(ma_score(probe_winds_indices{i_grp},:), 2), ...
        sqrt(nanvar(ma_score(probe_winds_indices{i_grp},:), [], 2)./...
        sum(~isnan(ma_score(probe_winds_indices{i_grp},:)),2)), ...
        'k');
    title(['Group: ', num2str(i_grp), ', Combo']);
    if ~self_scale
        axis(axis_ranges)
    end
end
f3 = gcf;
set(f3, 'Position', fig_pos);
union_avg_policy =...
    union_avg_policy(~isnan(union_avg_policy(:, 1)), :);
csvwrite([outpath, 'policy_union'], union_avg_policy);

diff_response_y = union_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(union_probe_diff_grp,1), 1), numel(union_probe_diff_grp), 1), ...
    nan(numel(union_probe_diff_grp), 1), ...
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])'];
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'policy_union_probe_diff_grp'], [diff_response_y, diff_design_X]);

figure; hold on;
errorbar([1 2 3 4], nanmean(union_probe_diff_grp), ...
    sqrt(nanvar(union_probe_diff_grp)./sum(~isnan(union_probe_diff_grp))), 'k.')
bar([1 2 3 4], nanmean(union_probe_diff_grp));
title('union')
axis([0 5 -.15 .1])
f3 = gcf;
set(f3, 'Position', bars_pos);

%% save policy score per trial for use in interaction test:

pre_inds = {26:75, 426:475, 826:875, 1826:1875};
probe_inds = {76:125, 476:525, 876:925, 1876:1925};

score_succ_by_trial_grp = cell(2,4); % {(pre-probe; probe) x group}
for i_grp = 1:size(score_succ_by_trial_grp,2)
     temp_pre = reshape(nanmean(...
         score_s{i_grp}(pre_inds{i_grp}, score_signal_window, :), 2),...
         size(score_s{i_grp}(pre_inds{i_grp}, :, :),1),...
         size(score_s{i_grp}(pre_inds{i_grp}, :, :),3));
     
     temp_probe = reshape(nanmean(...
         score_s{i_grp}(probe_inds{i_grp}, score_signal_window, :), 2),...
         size(score_s{i_grp}(probe_inds{i_grp}, :, :),1),...
         size(score_s{i_grp}(probe_inds{i_grp}, :, :),3));
     
     score_succ_by_trial_grp{1, i_grp} = temp_pre;
     score_succ_by_trial_grp{2, i_grp} = temp_probe;
end
% save('score_by_trial_grp.mat', 'score_succ_by_trial_grp');

score_fail_by_trial_grp = cell(2,4); % {(pre-probe; probe) x group}
for i_grp = 1:size(score_fail_by_trial_grp,2)
     temp_pre = reshape(nanmean(...
         score_f{i_grp}(pre_inds{i_grp}, score_signal_window, :), 2),...
         size(score_f{i_grp}(pre_inds{i_grp}, :, :),1),...
         size(score_f{i_grp}(pre_inds{i_grp}, :, :),3));
     
     temp_probe = reshape(nanmean(...
         score_f{i_grp}(probe_inds{i_grp}, score_signal_window, :), 2),...
         size(score_f{i_grp}(probe_inds{i_grp}, :, :),1),...
         size(score_f{i_grp}(probe_inds{i_grp}, :, :),3));
     
     score_fail_by_trial_grp{1, i_grp} = temp_pre;
     score_fail_by_trial_grp{2, i_grp} = temp_probe;
end
save('score_by_trial_grp.mat', 'score_succ_by_trial_grp', 'score_fail_by_trial_grp');

%% plot the deviation timecourse
plot_trajectory_deviation

%% plot the average deviation in response to the probe onset
plot_probe_response