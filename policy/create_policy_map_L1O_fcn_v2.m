function policy = create_policy_map_L1O_fcn_v2(tilt_dir, tilt_mag, traj_dir, traj_loc, len_range)
% function [s, u] = create_policy_map_L1O
%
% Create a policy map for every subject in every group by a leave-one-out
% ("L1O") method. Uses data from all subjects across all groups, minus the
% one subject for whom the policy is being computed, to creates a policy
% map. A policy map is a 3-D matrix of lateral bins, direction headings,
% and track locations. For each bin, heading, and location coordinate, a
% mean and variance of tilt direction from successful trials is saved. 
% This policy map is used to compare a single trajectory's tilt signal to 
% the distrubtion of titl signals from successful trajectories that had 
% matching states.
%
% David Huberdeau, 9/14/2018

% [tilt_dir, traj_dir, tilt_mag, traj_loc] = control_policy_v5_fcn(len_range);
% load('tilt_dir.mat')
% load('tilt_mag.mat')
% load('traj_dir.mat')
% load('traj_loc.mat')
load('exact_track_dist_full_v3.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat') % replaced with
% above, 9/27/17

% fall_off_succ_TH = .8;
fall_off_succ_TH = len_range(2) + .1;
sample_size_th = 10;

% grp_remove_sub = {[], [6], [2,3], [11 12 19]};
grp_remove_sub = {}; %9/14/18 - if subjects are to be removed - do it later.
% REMOVE_SUB = [2 4 11 19 20];

% REMOVE_SUB = [];

% loop_discard_TH = 200;
% temp_base = rad2deg(unwrap(deg2rad(temp_base)));
% traj_basis{i_day} = temp_base(:, sum(temp_base > loop_discard_TH, 1) < 1);

%%
% label successful trajectories all through segment
traj_succ = cell(size(traj_loc));

for k_grp = 1:length(traj_succ)
    traj_succ{k_grp} = nan(size(traj_loc{k_grp}));
    for i_sub = 1:size(traj_succ{k_grp},3)
        traj_succ{k_grp}(:,:,i_sub) = repmat(exact_track_dist_fall{k_grp}(:,i_sub)' > fall_off_succ_TH, size(traj_loc{k_grp},1), 1);
    end
end

%% compute Pr(visit location | succ)
% Use only best day trials:
day_inds = {...
    [], ...
    [], ...
    [], ...
    [1801:1875, 1926:2000]};

% Use all successfull trials:  commented out 11/21/2018:
% day_inds = {...
%     [1:75, 126:200], ...
%     [1:475, 526:1000], ...
%     [1:875, 926:1000], ...
%     [1:1875, 1926:2000]}; % exclude probe trials only. 


% day_inds = {[51:150], [251:350], [451:550], [651:750], [851:950],...
%     1000+[51:150], 1000+[251:350], 1000+[451:550], 1000+[651:750],
%     1000+[851:950]}; 

% grp_days = [1 5 5 10];
% p_loc_g_succ = cell(1,n_subs);

num_l_bins = 10;
L_BINS = linspace(0,1,num_l_bins+1);
num_d_bins = 20; %number of direction bins
D_BINS = linspace(-180, 180, num_d_bins + 1);
% num_locs = 20;
% num_locs = round(100*(len_range(2) - len_range(1)));
num_locs = ceil(100*(len_range(2) - len_range(1)));
% NOTE: TRACK LENGTHS USED ARE .25 TO .45

% for i_day = 1:grp_days(i_grp)

er = cell(1,1); i_er = 1;
policy = cell(1,4);
for i_grp = 1:4
    
    n_subs = size(tilt_dir{i_grp},3);
    u_loc = cell(2,n_subs);
    
    % moved * from here

    n = nan(num_l_bins,num_d_bins,num_locs);
    u = nan(num_l_bins,num_d_bins,num_locs);
    u_s = nan(num_l_bins,num_d_bins,num_locs);
    for i_sub = 1:n_subs
        succ_day = nan(100000, num_locs);
        loc_day = nan(100000, num_locs);
        tilt_day = nan(100000, num_locs);
        dir_day = nan(100000, num_locs);%num_locs was 20
        k = 1;
        for temp_grp = 1:length(day_inds)
            sub_inds = 1:size(traj_succ{temp_grp},3);
            % * moved to here
            if i_grp == temp_grp
                this_sub = i_sub;
            else
                this_sub = [];
            end
            succ_day_ = permute(traj_succ{temp_grp}(:,day_inds{temp_grp},setdiff(sub_inds, this_sub)), [2, 3, 1]);
            succ_day_t = reshape(succ_day_, size(succ_day_, 1)*size(succ_day_,2), size(succ_day_,3));
            succ_day(k:(k + size(succ_day_t,1) - 1), :) = succ_day_t;

            loc_day_ = permute(traj_loc{temp_grp}(:,day_inds{temp_grp},setdiff(sub_inds, this_sub)), [2, 3, 1]);
            loc_day_t = reshape(loc_day_, size(loc_day_, 1)*size(loc_day_,2), size(succ_day_,3));
            loc_day(k:(k + size(loc_day_t,1) - 1), :) = loc_day_t;

            tilt_day_ = permute(tilt_dir{temp_grp}(:,day_inds{temp_grp},setdiff(sub_inds, this_sub)), [2 3 1]);
            tilt_day_t = reshape(tilt_day_, size(tilt_day_, 1)*size(tilt_day_,2), size(succ_day_,3));
            tilt_day(k:(k + size(tilt_day_t,1) - 1), :) = tilt_day_t;

            dir_day_ = permute(traj_dir{temp_grp}(:,day_inds{temp_grp},setdiff(sub_inds, this_sub)), [2 3 1]);
            dir_day_t = reshape(dir_day_, size(dir_day_, 1)*size(dir_day_,2), size(succ_day_,3));
            dir_day(k:(k + size(dir_day_t,1) - 1), :) = dir_day_t;
            k = k + size(succ_day_t,1);
        end
        try
            temp_loc = loc_day(succ_day(:,1)>0,:);
            temp_tilt = tilt_day(succ_day(:,1)>0,:);
            temp_dir = dir_day(succ_day(:,1)>0,:);

            temp_n = nan(num_l_bins, num_d_bins, num_locs);
            for i_loc = 1:size(temp_loc,2)
                temp_n_ = hist3([temp_loc(:,i_loc), temp_dir(:, i_loc)], {L_BINS, D_BINS});
                temp_n(:,:,i_loc) = temp_n_(1:num_l_bins,1:num_d_bins);
            end
            n(:, :, :, i_sub) = temp_n;
            for i_loc = 1:num_locs
                for i_bin = 1:(num_l_bins-1)
                    for i_dir = 1:(num_d_bins-1)
                        crit1 = temp_loc(:,i_loc) > L_BINS(i_bin) & temp_loc(:,i_loc) <= L_BINS(i_bin + 1);
                        crit2 = temp_dir(:,i_loc) > D_BINS(i_dir) & temp_dir(:,i_loc) <= D_BINS(i_dir + 1);
                        if sum(crit1&crit2) >= sample_size_th
                            u(i_bin, i_dir, i_loc) = mean(temp_tilt(crit1 & crit2, i_loc));
                            u_s(i_bin, i_dir, i_loc) = var(temp_tilt(crit1 & crit2, i_loc));
                        end
                    end
                end
            end
        catch
            er{i_er} = lasterror;
            assignin('base', 'er', er);
            i_er = i_er + 1;
            warning('Error computing policy tilt and tilt variability.')
        end
        u_loc{1, i_sub} = u;
        u_loc{2, i_sub} = u_s;
    end
    policy{i_grp} = u_loc;
    disp(['Policy map complete']);
end