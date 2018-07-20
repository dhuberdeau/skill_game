% Plot the track with policy map collapsed across direction space and a few
% example trajectories.
%
% David Huberdeau 07/20/2018

%% load data
load('TracksStructure_v2.mat')
load('group_traj_x.mat')
load('group_traj_y.mat')
load('policy_all.mat')
load('deviation_maps_LEARN_SUCCESSES.mat')
load('deviation_maps_LEARN_FAILURES.mat')
load('exact_track_dist_full_v3.mat')

%% plot path
n_samp = size(p{4}{1,1}, 2); % samples is first dim
n_l = size(p{4}{1,1}, 1);  % lateral samples is second dim
n_sub = size(p{4}{1,1}, 3); % subjects is third dim

track_out = plotPath_fine(Tracks.fast_challenge_1.Paths); close;
dist_1 = (sqrt(diff(track_out(1, :)).^2 + diff(track_out(2, :)).^2));
dist_2 = (sqrt(diff(track_out(3, :)).^2 + diff(track_out(4, :)).^2));
cdist_1 = cumsum(dist_1)/sum(dist_1);
cdist_2 = cumsum(dist_2)/sum(dist_2);

% th_1 = .2; th_2 = .5; % th_2 = .46; %using 750 ms instead of end point
th_1 = .25; th_2 = .6;
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
plot([track_out(1,k_1_th); track_out(3,k_2_th)], [track_out(2,k_1_th); track_out(4,k_2_th)], 'r--')

for i_loc = 1:length(k_linsp_1)
    x1 = track_out(1, k_linsp_1(i_loc));
    x2 = track_out(3, k_linsp_2(i_loc));
    y1 = track_out(2, k_linsp_1(i_loc));
    y2 = track_out(4, k_linsp_2(i_loc));
    
    for i_l = 1:n_l
        s_x = x1 + (x2 - x1)/n_l*i_l;
        s_y = y1 + (y2 - y1)/n_l*i_l;

        plot(s_x, s_y, 'k.');
    end
end

%% Select success and failure trials

