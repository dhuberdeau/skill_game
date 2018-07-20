function [zu_trial_groups, zu_trial_probe_groups] = compute_policy_deviation_map(pre_wind, probe_wind, succ)
% function compute_policy_deviation_map
%
% Use the policy map (leave-one-out) to compute each subject's policy
% deviation map. This need only be done once, and additional analysis
% scripts can load the saved deviation maps.
%
% David Huberdeau 07/20/2018

%% load all data
load('tilt_dir.mat')
load('tilt_mag.mat')
load('traj_dir.mat')
load('traj_loc.mat')
load('exact_track_dist_full_v3.mat')
load('policy_all.mat');

%% define constants and parameters
fall_off_succ_TH = .8;
z_scr_TH = 50;

num_l_bins = 10;
L_BINS = linspace(0,1,num_l_bins+1);
num_d_bins = 20; %number of direction bins
D_BINS = linspace(-180, 180, num_d_bins + 1);
num_locs = 20;
grp_days = [1 5 5 10];

% pre_wind = {26:75, 226:275, 426:475, 626:675, 826:875, 1026:1075, 1226:1275, 1426:1475, 1626:1675, 1826:1875};
% probe_wind = {76:125, 276:325, 476:525, 676:725, 876:925, 1076:1125, 1276:1325, 1476:1525, 1676:1725, 1876:1925};
%% label successful trajectories all through segment
traj_succ = cell(size(traj_loc));

for k_grp = 1:length(traj_succ)
    traj_succ{k_grp} = nan(size(traj_loc{k_grp}));
    for i_sub = 1:size(traj_succ{k_grp},3)
        traj_succ{k_grp}(:,:,i_sub) = repmat(exact_track_dist_fall{k_grp}(:,i_sub)' > fall_off_succ_TH, size(traj_loc{k_grp},1), 1);
    end
end

%% define output variables
zu_trial_groups = cell(1,length(grp_days));
zu_trial_probe_groups = cell(1,length(grp_days));


%% compute Deviation map for pre-probe and probe windows

for i_grp = 1:4
    n_subs = size(exact_track_dist_fall{i_grp},2);
    u_loc = p{i_grp};
    zu_trial = nan(length(pre_wind{1}), num_locs, n_subs, grp_days(i_grp)); %10 days

    for i_day = 1:grp_days(i_grp)
        succ_day = permute(traj_succ{i_grp}(:,pre_wind{i_day},:), [2, 1, 3]);
        loc_day = permute(traj_loc{i_grp}(:,pre_wind{i_day},:), [2, 1, 3]);
        tilt_day = permute(tilt_dir{i_grp}(:,pre_wind{i_day},:), [2 1 3]);
        dir_day = permute(traj_dir{i_grp}(:,pre_wind{i_day},:), [2 1 3]);

        zu_day = nan(size(tilt_day));
    %     zus_day = nan(num_l_bins, num_d_bins, num_locs, n_subs);
        for i_sub = 1:size(succ_day,3)
            trial_inds = 1:size(succ_day,1);
            if succ
                % compute deviation map for ultimately successful trials
                target_trials = trial_inds(succ_day(:,1,i_sub)>0);
            else
                % compute deviation map for ultimately failure trials
                target_trials = trial_inds(succ_day(:,1,i_sub)<1);
            end

            zu_sub = nan(size(succ_day,1), size(succ_day,2));
            for i_tr = target_trials
                try
                    temp_loc_tr = loc_day(i_tr,:,i_sub);
                    temp_tilt_tr = tilt_day(i_tr,:,i_sub);
                    temp_dir_tr = dir_day(i_tr,:,i_sub);

                    for i_loc = 1:num_locs
                        for i_bin = 1:num_l_bins
                            for i_dir = 1:num_d_bins
                                crit1 = temp_loc_tr(i_loc) > L_BINS(i_bin) & temp_loc_tr(i_loc) <= L_BINS(i_bin + 1);
                                crit2 = temp_dir_tr(i_loc) > D_BINS(i_dir) & temp_dir_tr(i_loc) <= D_BINS(i_dir + 1);
        %                         temp_tilt_tr = temp_tilt(crit1 & crit2, i_loc);
        %                         zu_temp = nan(size(temp_tilt_tr,1), 1);
                                if crit1 && crit2
                                    mu = u_loc{1,i_sub}(i_bin, i_dir, i_loc);
    %                                 s_temp = reshape(u_conglom{2,i_sub}(i_bin, i_dir, i_loc, :), 10, 1); %10 for day
                                    sig = sqrt(u_loc{2,i_sub}(i_bin, i_dir, i_loc));
                                    if ~isnan(mu) && ~isnan(sig)
                                        if temp_tilt_tr(i_loc) > -180 && temp_tilt_tr(i_loc) < 180
                                            zscore_temp = abs((temp_tilt_tr(i_loc) - mu)./sig);
                                            if zscore_temp < z_scr_TH
        %                                     zscore_temp = ((temp_tilt_tr(i_loc) - mu));
        %                                     if zscore_temp < 180 & zscore_temp > -180
                                                zu_sub(i_tr, i_loc) = zscore_temp;
                                            else
                                                zu_sub(i_tr, i_loc) = nan;
                                            end
                                        else
                                            zu_sub(i_tr, i_loc) = nan;
                                        end
                                    else
                                        a=1+1;
                                    end
                                end
                            end
                        end
                    end
                catch er
                    warning('wtf');
                end
            end
        zu_day(:, :, i_sub) = zu_sub;    
        end
        zu_trial(:,:,:, i_day) = zu_day;
    %     z_u{1, i_day} = zu_day;
    %     z_u{2, i_day} = zus_day;
        disp(['Day: ', num2str(i_day)]);
    end
    zu_trial_groups{i_grp} = zu_trial;
    
    zu_trial_probe = nan(length(probe_wind{1}), num_locs, n_subs, grp_days(i_grp)); %10 days

    for i_day = 1:grp_days(i_grp)
        succ_day = permute(traj_succ{i_grp}(:,probe_wind{i_day},:), [2, 1, 3]);
        loc_day = permute(traj_loc{i_grp}(:,probe_wind{i_day},:), [2, 1, 3]);
        tilt_day = permute(tilt_dir{i_grp}(:,probe_wind{i_day},:), [2 1 3]);
        dir_day = permute(traj_dir{i_grp}(:,probe_wind{i_day},:), [2 1 3]);

        zu_day = nan(size(tilt_day));
    %     zus_day = nan(num_l_bins, num_d_bins, num_locs, n_subs);
        for i_sub = 1:size(succ_day,3)
            trial_inds = 1:size(succ_day,1);
            if succ
                % compute deviation map for ultimately successful trials
                target_trials = trial_inds(succ_day(:,1,i_sub)>0);
            else
                % compute deviation map for ultimately failure trials
                target_trials = trial_inds(succ_day(:,1,i_sub)<1);
            end

            zu_sub = nan(size(succ_day,1), size(succ_day,2));
            for i_tr = target_trials
                try
                    temp_loc_tr = loc_day(i_tr,:,i_sub);
                    temp_tilt_tr = tilt_day(i_tr,:,i_sub);
                    temp_dir_tr = dir_day(i_tr,:,i_sub);

                    for i_loc = 1:num_locs
                        for i_bin = 1:num_l_bins
                            for i_dir = 1:num_d_bins
                                crit1 = temp_loc_tr(:,i_loc) > L_BINS(i_bin) & temp_loc_tr(:,i_loc) <= L_BINS(i_bin + 1);
                                crit2 = temp_dir_tr(:,i_loc) > D_BINS(i_dir) & temp_dir_tr(:,i_loc) <= D_BINS(i_dir + 1);
                                if crit1 && crit2
                                    mu = u_loc{1,i_sub}(i_bin, i_dir, i_loc);
    %                                 s_temp = reshape(u_conglom{2,i_sub}(i_bin, i_dir, i_loc, :), 10, 1);
                                    sig = sqrt(u_loc{2,i_sub}(i_bin, i_dir, i_loc));
                                    if ~isnan(mu) && ~isnan(sig)
                                        if temp_tilt_tr(i_loc) > -180 && temp_tilt_tr(i_loc) < 180
                                            zscore_temp = abs((temp_tilt_tr(i_loc) - mu)./sig);
                                            if zscore_temp < z_scr_TH
        %                                     zscore_temp = ((temp_tilt_tr(i_loc) - mu));
        %                                     if zscore_temp < 180 & zscore_temp > -180
                                                zu_sub(i_tr, i_loc) = zscore_temp;
                                            else
                                                zu_sub(i_tr, i_loc) = nan;
                                            end
                                        else
                                            zu_sub(i_tr, i_loc) = nan;
                                        end
                                    else
                                       a = 1+1; 
                                    end
                                end
                            end
                        end
                    end
                catch
                    warning('wtf');
                end
            end
        zu_day(:, :, i_sub) = zu_sub;    
        end
        zu_trial_probe(:,:,:, i_day) = zu_day;
        disp(['Day: ', num2str(i_day)]);
    end
    zu_trial_probe_groups{i_grp} = zu_trial_probe;
end