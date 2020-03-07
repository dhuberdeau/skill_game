
%% define output data path
% outpath = ['/Users/david/OneDrive/Documents/JHU/BLAM_lab/',...
%     'Projects/skillLearning_antGame_A/Matlab/public_analysis/data/'];
outpath = pwd;

%% load data

load('traj_pos_25_90.mat');
load('exact_track_dist_full_v3.mat')
load('subject_basis.mat');

%% rescale distance travelled to track lengths (not arbitrary game units)
track_len_conv_fct = 226.8773;
for i_grp = 1:4
    traj_pos{i_grp} = traj_pos{i_grp}./track_len_conv_fct;
end

%% define some parameters:
n_samps = 45; % good balance of trajectory length and available samples
ind_samps = [1:n_samps, 90 + (1:n_samps)]; % [x, y] position
wind_size = 25; %trials in a window.
non_nan_th = 7; %min. samples from 1 subject to include them.

%% cancel subjects if needed:
sub_grp = [];
sub_num = []; %better to exclude subjects at last stage before stats, not here

for i_sub = 1:length(sub_num)
    traj_pos{sub_grp(i_sub)}(:,:,sub_num(i_sub)) = nan;
end

%% get scores across windows of trials and the difference between pre-probe and probe.
k_pc = 1;
score_mats = {nan(200, 21), nan(1000, 20), nan(1000, 20), nan(2000, 20)};
for i_grp = 1:4
    for i_sub = 1:size(traj_pos{i_grp}, 3)
        data_mat = traj_pos{i_grp}(ind_samps, :, i_sub);
        data_mat_ = data_mat - repmat(subject_basis_mean{i_grp, i_sub}, 1, size(data_mat,2));
        score = subject_basis{i_grp, i_sub}(:, k_pc)'*data_mat_;
        score_mats{i_grp}(:, i_sub) = score';
    end
end

% plot variances:
var_mats = {nan(200/wind_size, 21), nan(1000/wind_size, 20), nan(1000/wind_size, 20), nan(2000/wind_size,20)};
men_mats = {nan(200/wind_size, 21), nan(1000/wind_size, 20), nan(1000/wind_size, 20), nan(2000/wind_size,20)};
abs_men_mats = {nan(200/wind_size, 21), nan(1000/wind_size, 20), nan(1000/wind_size, 20), nan(2000/wind_size,20)};
for i_grp = 1:4
    temp_dat = score_mats{i_grp};
    varss = nan(size(temp_dat,1)/wind_size, size(temp_dat,2));
    menss = nan(size(temp_dat,1)/wind_size, size(temp_dat,2));
    abs_menss = nan(size(temp_dat,1)/wind_size, size(temp_dat,2));
    for i_sub = 1:size(temp_dat,2)
        sub_scores = temp_dat(:, i_sub);
        reshp_scores = reshape(sub_scores, wind_size, size(sub_scores,1)/wind_size);

        valid_winds = sum(~isnan(reshp_scores),1) >= non_nan_th;

        % fill in valid variances:
        sub_vars = nan(1, length(valid_winds));
        temp_vars_ = nanvar(reshp_scores);
        sub_vars(valid_winds) = temp_vars_(valid_winds);
        varss(:, i_sub) = sub_vars';

        % fill in valid means:
        sub_mens = nan(1, length(valid_winds));
        temp_mens_ = nanmean(reshp_scores);
        sub_mens(valid_winds) = temp_mens_(valid_winds);
%         menss(:, i_sub) = abs(sub_mens)';
        menss(:, i_sub) = (sub_mens)';
        
        % fill in valid absolute-value means:
        sub_mens = nan(1, length(valid_winds));
        temp_mens_ = nanmean(abs(reshp_scores));
        sub_mens(valid_winds) = temp_mens_(valid_winds);
%         menss(:, i_sub) = abs(sub_mens)';
        abs_menss(:, i_sub) = (sub_mens)';
    end
    var_mats{i_grp} = varss;
    men_mats{i_grp} = menss;
    abs_men_mats{i_grp} = abs_menss;
end

%% plot means and variances over windows:
training_winds_grp = {[1 1 1 0 0 1 1 1],...
    [ones(1, 19), 0, 0, ones(1, 19)],...
    [ones(1, 35), 0, 0, ones(1, 3)],...
    [ones(1, 40 + 35), 0, 0, ones(1, 3)]};
clrs = {'y', 'r', 'b', 'g'};

% plot vars:
figure;
for i_grp = 1:4
    subplot(1,4,i_grp); hold on;
    inds = 1:size(var_mats{i_grp},1);
    errorfield(inds, nanmean(sqrt(var_mats{i_grp}),2), ...
        sqrt(nanvar(var_mats{i_grp}, [], 2)./size(var_mats{i_grp},2)), clrs{i_grp});
    errorfield(inds(~training_winds_grp{i_grp}), nanmean(sqrt(var_mats{i_grp}(~training_winds_grp{i_grp}, :)),2), ...
        sqrt(nanvar(var_mats{i_grp}(~training_winds_grp{i_grp}, :), [], 2)./size(var_mats{i_grp}(~training_winds_grp{i_grp}, :),2)), 'k');
    axis([0 80 .045 .08])
end


% plot means:
figure;
for i_grp = 1:4
    subplot(1,4,i_grp); hold on;
    inds = 1:size(men_mats{i_grp},1);
    errorfield(inds, nanmean((men_mats{i_grp}),2), ...
        sqrt(nanvar(men_mats{i_grp}, [], 2)./size(men_mats{i_grp},2)), clrs{i_grp});
    errorfield(inds(~training_winds_grp{i_grp}), nanmean((men_mats{i_grp}(~training_winds_grp{i_grp}, :)),2), ...
        sqrt(nanvar(men_mats{i_grp}(~training_winds_grp{i_grp}, :), [], 2)./size(men_mats{i_grp}(~training_winds_grp{i_grp}, :),2)), 'k');
    axis([0 80 -.05 .05])
end

%% save data for analysis in R
training_winds_grp = {boolean([1 1 1 0 0 1 1 1]),...
    boolean([ones(1, 19), 0, 0, ones(1, 19)]),...
    boolean([ones(1, 35), 0, 0, ones(1, 3)]),...
    boolean([ones(1, 40 + 35), 0, 0, ones(1, 3)])};
response_y = nan(4000, 1); % pre-alocate matrix with sufficient size
predictor_X = nan(4000, 3); %use window and training (binar) as predictors
k_samp = 1; k_abs_sub = 1; %each subject is different, dont use i_sub index
for i_grp = 1:length(var_mats)
    for i_sub = 1:size(var_mats{i_grp}, 2)
        resp_training = var_mats{i_grp}(training_winds_grp{i_grp}, i_sub);
        wind_inds = 1:size(var_mats{i_grp},1);
        
        response_y(k_samp - 1 + (1:length(resp_training))) = resp_training;
        predictor_X(k_samp - 1 + (1:length(resp_training)), 1) = ...
            wind_size*(wind_inds(training_winds_grp{i_grp}));
        predictor_X(k_samp - 1 + (1:length(resp_training)), 2) = ...
            i_grp*ones(length(resp_training),1);
        predictor_X(k_samp - 1 + (1:length(resp_training)), 3) = ...
            k_abs_sub*ones(length(resp_training), 1);
        k_samp = k_samp + length(resp_training);
        k_abs_sub = k_abs_sub + 1;
    end
end
predictor_X = predictor_X(~isnan(response_y), :);
response_y = response_y(~isnan(response_y));
csvwrite([outpath, 'variability_data_pc', num2str(k_pc)], [response_y, predictor_X]);

response_y = nan(4000, 1); % pre-alocate matrix with sufficient size
predictor_X = nan(4000, 3); %use window and training (binar) as predictors
k_samp = 1; k_abs_sub = 1;
for i_grp = 1:length(men_mats)
    for i_sub = 1:size(men_mats{i_grp}, 2)
        resp_training = men_mats{i_grp}(training_winds_grp{i_grp}, i_sub);
        wind_inds = 1:size(men_mats{i_grp},1);
        
        response_y(k_samp - 1 + (1:length(resp_training))) = resp_training;
        predictor_X(k_samp - 1 + (1:length(resp_training)), 1) = ...
            wind_size*(wind_inds(training_winds_grp{i_grp}));
        predictor_X(k_samp - 1 + (1:length(resp_training)), 2) = ...
            i_grp*ones(length(resp_training),1);
        predictor_X(k_samp - 1 + (1:length(resp_training)), 3) = ...
            k_abs_sub*ones(length(resp_training), 1);
        k_samp = k_samp + length(resp_training);
        k_abs_sub = k_abs_sub + 1;
    end
end
predictor_X = predictor_X(~isnan(response_y), :);
response_y = response_y(~isnan(response_y));
csvwrite([outpath, 'mean_data_pc', num2str(k_pc)], [response_y, predictor_X]);


response_y = nan(4000, 1); % pre-alocate matrix with sufficient size
predictor_X = nan(4000, 3); %use window and training (binar) as predictors
k_samp = 1; k_abs_sub = 1;
for i_grp = 1:length(abs_men_mats)
    for i_sub = 1:size(abs_men_mats{i_grp}, 2)
        resp_training = abs_men_mats{i_grp}(training_winds_grp{i_grp}, i_sub);
        wind_inds = 1:size(abs_men_mats{i_grp},1);
        
        response_y(k_samp - 1 + (1:length(resp_training))) = resp_training;
        predictor_X(k_samp - 1 + (1:length(resp_training)), 1) = ...
            wind_size*(wind_inds(training_winds_grp{i_grp}));
        predictor_X(k_samp - 1 + (1:length(resp_training)), 2) = ...
            i_grp*ones(length(resp_training),1);
        predictor_X(k_samp - 1 + (1:length(resp_training)), 3) = ...
            k_abs_sub*ones(length(resp_training), 1);
        k_samp = k_samp + length(resp_training);
        k_abs_sub = k_abs_sub + 1;
    end
end
predictor_X = predictor_X(~isnan(response_y), :);
response_y = response_y(~isnan(response_y));
csvwrite([outpath, 'abs_mean_data_pc', num2str(k_pc)], [response_y, predictor_X]);
%%
% prep mean data for export to R for statistical analyses
mean_probe_diff_grp = nan(21, 4);
for i_grp = 1:4
    inds = 1:size(abs_men_mats{i_grp},1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;

    probe_winds = nanmean(abs_men_mats{i_grp}(probe_inds, :));
    pre_p_winds = nanmean(abs_men_mats{i_grp}(pre_p_inds, :));

    mean_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
end

diff_response_y = mean_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(mean_probe_diff_grp,1), 1), numel(mean_probe_diff_grp), 1), ... %groups
    k_pc*ones(numel(mean_probe_diff_grp), 1), ... % PC
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])']; % subjects
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'abs_mean_probe_difference_pc', num2str(k_pc)], [diff_response_y, diff_design_X]);

%%
% prep mean data for export to R for statistical analyses
mean_probe_diff_grp = nan(21, 4);
for i_grp = 1:4
    inds = 1:size(men_mats{i_grp},1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;

    probe_winds = nanmean(men_mats{i_grp}(probe_inds, :));
    pre_p_winds = nanmean(men_mats{i_grp}(pre_p_inds, :));

    mean_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
end

diff_response_y = mean_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(mean_probe_diff_grp,1), 1), numel(mean_probe_diff_grp), 1), ... %groups
    k_pc*ones(numel(mean_probe_diff_grp), 1), ... % PC
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])']; % subjects
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'mean_probe_difference_pc', num2str(k_pc)], [diff_response_y, diff_design_X]);

figure; hold on;
errorbar([1 2 3 4], nanmean(mean_probe_diff_grp), ...
    sqrt(nanvar(mean_probe_diff_grp)./sum(~isnan(mean_probe_diff_grp))), 'k.')
bar([1 2 3 4], nanmean(mean_probe_diff_grp));

%%
% prep variability data for export to R for statistical analyses
var_probe_diff_grp = nan(21, 4);
for i_grp = 1:4
    inds = 1:size(var_mats{i_grp},1);
    probe_inds = inds(~training_winds_grp{i_grp});
    pre_p_inds = probe_inds - 2;

    probe_winds = nanmean(var_mats{i_grp}(probe_inds, :));
    pre_p_winds = nanmean(var_mats{i_grp}(pre_p_inds, :));

    var_probe_diff_grp(1:length(probe_winds), i_grp) = (probe_winds - pre_p_winds)';
end

diff_response_y = var_probe_diff_grp(:);
diff_design_X = [...
    reshape(repmat((1:4), size(var_probe_diff_grp,1), 1), numel(var_probe_diff_grp), 1), ...
    k_pc*ones(numel(var_probe_diff_grp), 1), ...
    ([1:21, 21+(1:20), nan, 41+(1:20), nan, 61+(1:20), nan])'];
response_nan_inds = ~isnan(diff_response_y);
diff_response_y = diff_response_y(response_nan_inds);
diff_design_X = diff_design_X(response_nan_inds, :);
csvwrite([outpath, 'var_probe_difference_pc', num2str(k_pc)], [diff_response_y, diff_design_X]);

figure; hold on;
errorbar([1 2 3 4], nanmean(var_probe_diff_grp), ...
    sqrt(nanvar(var_probe_diff_grp)./sum(~isnan(var_probe_diff_grp))), 'k.')
bar([1 2 3 4], nanmean(var_probe_diff_grp));