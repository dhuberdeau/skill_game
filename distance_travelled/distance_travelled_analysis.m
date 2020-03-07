
%% define output data path
outpath = ['/Users/david/OneDrive/Documents/JHU/BLAM_lab/',...
    'Projects/skillLearning_antGame_A/Matlab/public_analysis/data/'];

%% load data
load('exact_track_dist_full_v3.mat')

%% define some parameters:
wind_size = 25; %trials in a window.
axis_dims = [0 2000 0.25 0.8];
%% cancel subjects if needed:
sub_grp = [];
sub_num = []; %better to exclude subjects at last stage before stats, not here

for i_sub = 1:length(sub_num)
    traj_pos{sub_grp(i_sub)}(:,:,sub_num(i_sub)) = nan;
end

%% get scores across windows of trials and the difference between pre-probe and probe.
training_winds_grp = {boolean([1 1 1 0 0 1 1 1]),...
    boolean([ones(1, 19), 0, 0, ones(1, 19)]),...
    boolean([ones(1, 35), 0, 0, ones(1, 3)]),...
    boolean([ones(1, 40 + 35), 0, 0, ones(1, 3)])};
clrs = {'y', 'r', 'b', 'g'};
windows_cell = ...
    {nan(length(training_winds_grp{1}), size(exact_track_dist_fall{1}, 1)),...
    nan(length(training_winds_grp{2}), size(exact_track_dist_fall{2}, 1)),...
    nan(length(training_winds_grp{3}), size(exact_track_dist_fall{3}, 1)),...
    nan(length(training_winds_grp{4}), size(exact_track_dist_fall{4}, 1))};

difference_cell = ...
    {nan(length(training_winds_grp{1}), size(exact_track_dist_fall{1}, 1)),...
    nan(length(training_winds_grp{2}), size(exact_track_dist_fall{2}, 1)),...
    nan(length(training_winds_grp{3}), size(exact_track_dist_fall{3}, 1)),...
    nan(length(training_winds_grp{4}), size(exact_track_dist_fall{4}, 1))};

data_window = nan(50000, 4); % total samples x variables
data_diff_probe = nan(50000, 4);
k_wind = 1;
k_prob = 1;
k_sub = 0;

figure; 
for i_grp = 1:4
    subplot(1,4,i_grp); hold on;

    group_dist = exact_track_dist_fall{i_grp};

    dist_ma = mafilt(group_dist, wind_size, 1);
        
    dist_inds = wind_size*(1:size(dist_ma,1)); 
    dist_mean = nanmean(dist_ma, 2);
    dist_se = sqrt(nanvar(dist_ma, [], 2)./sum(~isnan(dist_ma),2));
    errorfield(dist_inds, dist_mean, dist_se, clrs{i_grp});
    
    errorfield(...
        dist_inds(~training_winds_grp{i_grp}),...
        dist_mean(~training_winds_grp{i_grp}),...
        dist_se(~training_winds_grp{i_grp}),...
        'k-');
    
    axis(axis_dims);

    % collect data for window analysis:
    dist_train = dist_ma(training_winds_grp{i_grp}, :);
    dist_wind = repmat(dist_inds(training_winds_grp{i_grp})', 1, size(dist_ma,2));
    subjs = repmat(k_sub + (1:size(dist_ma,2)), size(dist_train,1), 1);
    data_window(k_wind:(k_wind - 1 + numel(dist_train)), :) = ...
        [dist_train(:), dist_wind(:), i_grp*ones(numel(dist_train), 1), subjs(:)];
    k_wind = k_wind + numel(dist_train);
    
    % collect data for probe difference analysis:
    inds = 1:size(dist_ma,1);
    inds_probe = inds(~training_winds_grp{i_grp});
    inds_preprobe = inds_probe - 2;
    dist_probe = nanmean(dist_ma(inds_probe, :),1);
    dist_preprobe = nanmean(dist_ma(inds_preprobe, :),1);
    data_diff_probe(k_prob:(k_prob - 1 + length(dist_probe)), :) = ...
        [(dist_probe - dist_preprobe)', ...
        i_grp*ones(length(dist_probe), 1),...
        nan(length(dist_probe), 1),...
        k_sub + (1:length(dist_probe))'];
    k_prob = k_prob + length(dist_probe);
    
    k_sub = k_sub + size(dist_ma,2);
end
valid_samps = ~isnan(data_window(:, 1));
data_window = data_window(valid_samps, :);
csvwrite([outpath, 'dist_travelled_window'], data_window);

valid_samps = ~isnan(data_diff_probe(:,1));
data_diff_probe = data_diff_probe(valid_samps, :);
csvwrite([outpath, 'dist_travelled_probe_diff'], data_diff_probe);