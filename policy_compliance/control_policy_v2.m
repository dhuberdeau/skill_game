%% load data
% load('exact_track_dist_full_strictExclusions_v2.mat')
load('exact_track_dist_full_v3.mat')
load('TracksStructure_v2.mat')
load('group_traj_x.mat')
load('group_traj_y.mat')
% load('group_advanced_fall_data.mat')
load('falloff_traj_ind.mat')
load('group_tilt_dir.mat')
load('group_tilt_mag.mat')
load('group_track_num.mat') % 1 = _1, 0 = _mirror

%% 

n_samp = 10;%20

track_out = plotPath_fine(Tracks.fast_challenge_1.Paths); close;
dist_1 = (sqrt(diff(track_out(1, :)).^2 + diff(track_out(2, :)).^2));
dist_2 = (sqrt(diff(track_out(3, :)).^2 + diff(track_out(4, :)).^2));
cdist_1 = cumsum(dist_1)/sum(dist_1);
cdist_2 = cumsum(dist_2)/sum(dist_2);

% th_1 = .2; th_2 = .5; % th_2 = .46; %using 750 ms instead of end point
th_1 = .25; th_2 = .45;
[junk, k_min_30_1] = min(abs(cdist_1 - th_1));
[junk, k_min_30_2] = min(abs(cdist_2 - th_1));

linsp = linspace(th_1, th_2, n_samp);
k_linsp_1 = nan(1,n_samp);
k_linsp_2 = nan(1,n_samp);
for i_l = 1:length(linsp)
    [junk, k_linsp_1(i_l)] = min(abs(linsp(i_l) - cdist_1));
end
for i_l = 1:length(linsp)
    [junk, k_linsp_2(i_l)] = min(abs(linsp(i_l) - cdist_2));
end

[junk, k_1_th] = min(abs(th_1 - cdist_1));
[junk, k_2_th] = min(abs(th_1 - cdist_2));

figure
hold on;
plot(track_out(1,:), track_out(2,:))
plot(track_out(3,:), track_out(4,:))
axis equal
plot([track_out(1,k_linsp_1); track_out(3,k_linsp_2)], [track_out(2,k_linsp_1); track_out(4,k_linsp_2)], 'k-');
plot([track_out(1,k_1_th); track_out(3,k_2_th)], [track_out(2,k_1_th); track_out(4,k_2_th)], 'r-')


%%
% take position (distance from left side), direction heading, tilt
% direction and tilt magnitude at
% each of equally-spaced rungs along track.

tilt_mag = {nan(n_samp, 200, 21), nan(n_samp, 1000, 20), nan(n_samp, 1000, 20), nan(n_samp, 2000, 20)};
tilt_dir = {nan(n_samp, 200, 21), nan(n_samp, 1000, 20), nan(n_samp, 1000, 20), nan(n_samp, 2000, 20)};
traj_dir = {nan(n_samp, 200, 21), nan(n_samp, 1000, 20), nan(n_samp, 1000, 20), nan(n_samp, 2000, 20)};
traj_loc = {nan(n_samp, 200, 21), nan(n_samp, 1000, 20), nan(n_samp, 1000, 20), nan(n_samp, 2000, 20)};
% trial_dt = {nan(200, 21), nan(n_samp, 1000, 20), nan(n_samp, 1000, 20), nan(n_samp, 2000, 20)}; %distance travelled (even beyond analysis window)

for i_grp = 1:4
    for i_sub = 1:size(group_traj_x{i_grp},3)
        for i_tr = 1:size(group_traj_x{i_grp},2)
            
            keep_inds = 1:falloff_traj_ind{i_grp}(i_tr, i_sub);
            
            traj_x = group_traj_x{i_grp}(keep_inds, i_tr, i_sub);
            traj_y = group_traj_y{i_grp}(keep_inds, i_tr, i_sub);
            tilt_d = rad2deg(unwrap(deg2rad(group_tilt_dir{i_grp}(keep_inds, i_tr, i_sub))));
            tilt_m = group_tilt_mag{i_grp}(keep_inds, i_tr, i_sub);
%             dist_t = exact_track_dist_fall{i_grp}(i_tr, i_sub);
            
            if group_track_num{i_grp}(i_tr, i_sub) == 0
                %mirror track
                traj_x = -traj_x;
                tilt_d = -tilt_d;
            end
            
            for i_samp = 1:length(k_linsp_1)
                traj_inds = 1:length(traj_x);
                int_out = lineSegmentIntersect(...
                    [track_out(1, k_linsp_1(i_samp)), track_out(2, k_linsp_1(i_samp)), track_out(3, k_linsp_2(i_samp)), track_out(4, k_linsp_2(i_samp))], ...
                    [traj_x(1:(end-1)), traj_y(1:(end-1)), traj_x(2:end), traj_y(2:end)]);
                intx_ind = traj_inds(int_out.intAdjacencyMatrix > 0);
                try
                    intx_ind = intx_ind(1);

                    traj_loc{i_grp}(i_samp, i_tr, i_sub) = int_out.intNormalizedDistance1To2(1, intx_ind);%sqrt((traj_y(intx_ind) - ).^2 + ().^2);
                    
                    traj_dir_all = atan2d(diff(traj_y((intx_ind - 1):intx_ind)), diff(traj_x((intx_ind - 1):intx_ind))) - 90;
                    traj_dir{i_grp}(i_samp, i_tr, i_sub) = traj_dir_all;

                    % No filtering of trajectories needed

                    tilt_dir_filt = sgolayfilt(tilt_d, 3, 5);

                    tilt_dir{i_grp}(i_samp, i_tr, i_sub) = tilt_dir_filt(intx_ind);
                    
                    tilt_mag_filt = sgolayfilt(tilt_m, 3, 5);

                    tilt_mag{i_grp}(i_samp, i_tr, i_sub) = tilt_mag_filt(intx_ind);
                    
                catch
                    intx_ind = nan;
                end
            end
        end
    end
end

%% save tilt_dir

save tilt_dir tilt_dir
save traj_dir traj_dir
save tilt_mag tilt_mag
save traj_loc traj_loc