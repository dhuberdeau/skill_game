function zu_trial_groups = compute_policy_deviation_map_fcn(policy, exact_track_dist_fall, len_range, succ, tilt_dir, traj_dir, tilt_mag, traj_loc)
% function compute_policy_deviation_map
%
% Use the policy map (leave-one-out) to compute each subject's policy
% deviation map. This need only be done once, and additional analysis
% scripts can load the saved deviation maps.
%
% David Huberdeau 09/14/2018

%% define parameters:
fall_off_succ_TH = len_range(2) + .1;

%% load all data

% f1=exist(['tilt_dir', num2str(round(10*len_range(end))), '.mat'], 'file');
% f2=exist(['traj_dir', num2str(round(10*len_range(end))), '.mat'], 'file');
% f3=exist(['tilt_mag', num2str(round(10*len_range(end))), '.mat'], 'file');
% f4=exist(['traj_loc', num2str(round(10*len_range(end))), '.mat'], 'file');
% if f1==2 && f2==2 && f3==2 && f4==2
%     load(['tilt_dir', num2str(round(10*len_range(end))), '.mat'])
%     load(['traj_dir', num2str(round(10*len_range(end))), '.mat'])
%     load(['tilt_mag', num2str(round(10*len_range(end))), '.mat'])
%     load(['traj_loc', num2str(round(10*len_range(end))), '.mat'])
% else
%     [tilt_dir, traj_dir, tilt_mag, traj_loc] = control_policy_v5_fcn(len_range);
%     save(['tilt_dir', num2str(round(10*len_range(end)))], 'tilt_dir');
%     save(['traj_dir', num2str(round(10*len_range(end)))], 'traj_dir');
%     save(['tilt_mag', num2str(round(10*len_range(end)))], 'tilt_mag');
%     save(['traj_loc', num2str(round(10*len_range(end)))], 'traj_loc');
% end


%% define constants and parameters

z_scr_TH = 10;

num_l_bins = 10;
L_BINS = linspace(0,1,num_l_bins+1);
num_d_bins = 20; %number of direction bins
D_BINS = linspace(-180, 180, num_d_bins + 1);
num_locs = round(100*(len_range(2) - len_range(1)));
grp_days = [1 5 5 10];

%% label successful trajectories all through segment
traj_succ = cell(size(traj_loc));

for k_grp = 1:length(traj_succ)
    traj_succ{k_grp} = nan(size(traj_loc{k_grp}));
    for i_sub = 1:size(traj_succ{k_grp},3)
        traj_succ{k_grp}(:,:,i_sub) = repmat(exact_track_dist_fall{k_grp}(:,i_sub)' > fall_off_succ_TH, size(traj_loc{k_grp},1), 1);
    end
end

%% define output variables
zu_trial_groups = cell(1,length(tilt_dir));

%% compute Deviation map for pre-probe and probe windows

for i_grp = 1:4
    u_loc = policy{i_grp};
    succ_grp = permute(traj_succ{i_grp}, [2, 1, 3]);
    loc_grp = permute(traj_loc{i_grp}, [2, 1, 3]);
    tilt_grp = permute(tilt_dir{i_grp}, [2 1 3]);
    dir_grp = permute(traj_dir{i_grp}, [2 1 3]);

    zu_trial = nan(size(tilt_grp));
    for i_sub = 1:size(succ_grp,3)
        trial_inds = 1:size(succ_grp,1);
        if succ == 1
            % compute deviation map for ultimately successful trials
            target_trials = trial_inds(succ_grp(:,1,i_sub)>0);
        elseif succ == 0
            % compute deviation map for ultimately failure trials
            target_trials = trial_inds(succ_grp(:,1,i_sub)<1);
        elseif succ == 2 
            % both failures and successes together
            target_trials = trial_inds;
        else
            error('Must specify subset of trials to include in analysis')
        end

        zu_sub = nan(size(succ_grp,1), size(succ_grp,2));
        for i_tr = target_trials
            try
                temp_loc_tr = loc_grp(i_tr,:,i_sub);
                temp_tilt_tr = tilt_grp(i_tr,:,i_sub);
                temp_dir_tr = dir_grp(i_tr,:,i_sub);

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
            catch err_pol_dev
                warning(err_pol_dev.message);
                warning([num2str(err_pol_dev.stack(1).line), ' : ', err_pol_dev.stack(1).file]);
            end
        end
        zu_trial(:, :, i_sub) = zu_sub;    
    end
    zu_trial_groups{i_grp} = zu_trial;
end