
%% define output data path
% outpath = ['/Users/david/OneDrive/Documents/JHU/BLAM_lab/',...
%     'Projects/skillLearning_antGame_A/Matlab/public_analysis/data/'];
outpath = pwd;

%% load data
load('traj_pos_25_90.mat');
load('exact_track_dist_full_v3.mat')

track_len_conv_fct = 226.8773;
for i_grp = 1:4
    traj_pos{i_grp} = traj_pos{i_grp}./track_len_conv_fct;
end

%% 
% collect all participants except the one of interest to use as "training" 
% data, and do a PCA using all of that data.
n_samps = 45;
ind_samps = [1:n_samps, 90 + (1:n_samps)];
subject_basis = cell(4, 21);
subject_basis_mean = cell(4,21);
subject_singularv = cell(4, 21);
grp_subs = {1:21, 1:20, 1:20, [1:11, 13:20]};
train_mat = nan(length(ind_samps), 10000);
for i_grp = 1:4
    for i_sub = 1:size(traj_pos{i_grp}, 3)
        k_samp = 1;
        for k_grp = 1:4
            for k_sub = 1:size(traj_pos{k_grp}, 3)
                this_sub_pos = traj_pos{k_grp}(ind_samps, :, k_sub);
                this_sub_pos_succ = this_sub_pos(:, ~isnan(this_sub_pos(n_samps, :)));
                train_mat(:, k_samp - 1 + (1:size(this_sub_pos_succ,2))) = ...
                    this_sub_pos_succ;
                k_samp = k_samp + size(this_sub_pos_succ,2);
            end
        end
        train_mat_0 = nanmean(train_mat, 2);
        train_mat_ = train_mat - repmat(train_mat_0, 1, size(train_mat,2));
        [u,s,v] = svd(train_mat_, 'econ');
        subject_basis{i_grp, i_sub} = u;
        subject_basis_mean{i_grp, i_sub} = train_mat_0;
        subject_singularv{i_grp, i_sub} = diag(s);
    end
end
save([outpath, 'subject_basis'],...
    'subject_basis', 'subject_basis_mean', 'subject_singularv');