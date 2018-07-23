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
load('falloff_traj_ind.mat')
load('group_tilt_dir.mat')

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

p__ = p{4}{1,10};
p_ = reshape(nanmean(p__, 2), 10, 20);

%% Failures

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
        s_x = x1 + (x2 - x1)/n_l*i_l - .5*(x2 - x1)/n_l;
        s_y = y1 + (y2 - y1)/n_l*i_l - .5*(y2 - y1)/n_l;

        plot(s_x, s_y, 'k.');
        
        
        plot(s_x + [0 2*cosd(90+p_(i_l, i_loc))], s_y + [0 2*sind(90+p_(i_l, i_loc))], 'r-')
    end
end

% Select success and failure trials
k_sub = 11;
inds = 1:2000;
fs = exact_track_dist_fall{4}(:, k_sub) < .8 & exact_track_dist_fall{4}(:, k_sub) > .45;
fsi = inds(fs);

k_samp = 1;
trajx = group_traj_x{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
trajy = group_traj_y{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
tilt = group_tilt_dir{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
plot(trajx, trajy, 'm-', 'LineWidth', 1)
for i_loc = 1:length(k_linsp_1)
    x1 = track_out(1, k_linsp_1(i_loc));
    x2 = track_out(3, k_linsp_2(i_loc));
    y1 = track_out(2, k_linsp_1(i_loc));
    y2 = track_out(4, k_linsp_2(i_loc));
    XY1 = [x1 y1 x2 y2];
    XY2 = [trajx(1:(end-1)), trajy(1:(end-1)), trajx(2:end), trajy(2:end)];
    intx = lineSegmentIntersect(XY1, XY2);
    plot(intx.intMatrixX(intx.intAdjacencyMatrix) + [0 2*cosd(90+tilt(intx.intAdjacencyMatrix))],...
        intx.intMatrixY(intx.intAdjacencyMatrix) + [0 2*sind(90+tilt(intx.intAdjacencyMatrix))],...
        'r-', 'LineWidth', 2)
end

k_samp = 2;
trajx = group_traj_x{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
trajy = group_traj_y{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
tilt = group_tilt_dir{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
plot(trajx, trajy, 'm-', 'LineWidth', 1)
for i_loc = 1:length(k_linsp_1)
    x1 = track_out(1, k_linsp_1(i_loc));
    x2 = track_out(3, k_linsp_2(i_loc));
    y1 = track_out(2, k_linsp_1(i_loc));
    y2 = track_out(4, k_linsp_2(i_loc));
    XY1 = [x1 y1 x2 y2];
    XY2 = [trajx(1:(end-1)), trajy(1:(end-1)), trajx(2:end), trajy(2:end)];
    intx = lineSegmentIntersect(XY1, XY2);
    plot(intx.intMatrixX(intx.intAdjacencyMatrix) + [0 2*cosd(90+tilt(intx.intAdjacencyMatrix))],...
        intx.intMatrixY(intx.intAdjacencyMatrix) + [0 2*sind(90+tilt(intx.intAdjacencyMatrix))],...
        'r-', 'LineWidth', 2)
end

%% Successes

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
        s_x = x1 + (x2 - x1)/n_l*i_l - .5*(x2 - x1)/n_l;
        s_y = y1 + (y2 - y1)/n_l*i_l - .5*(y2 - y1)/n_l;

        plot(s_x, s_y, 'k.');
        
        
        plot(s_x + [0 2*cosd(90+p_(i_l, i_loc))], s_y + [0 2*sind(90+p_(i_l, i_loc))], 'r-')
    end
end

fs = exact_track_dist_fall{4}(:, k_sub) > .6;
fsi = inds(fs);

k_samp = 1;
trajx = group_traj_x{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
trajy = group_traj_y{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
tilt = group_tilt_dir{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
plot(trajx, trajy, 'm-', 'LineWidth', 1)
for i_loc = 1:length(k_linsp_1)
    x1 = track_out(1, k_linsp_1(i_loc));
    x2 = track_out(3, k_linsp_2(i_loc));
    y1 = track_out(2, k_linsp_1(i_loc));
    y2 = track_out(4, k_linsp_2(i_loc));
    XY1 = [x1 y1 x2 y2];
    XY2 = [trajx(1:(end-1)), trajy(1:(end-1)), trajx(2:end), trajy(2:end)];
    intx = lineSegmentIntersect(XY1, XY2);
    plot(intx.intMatrixX(intx.intAdjacencyMatrix) + [0 2*cosd(90+tilt(intx.intAdjacencyMatrix))],...
        intx.intMatrixY(intx.intAdjacencyMatrix) + [0 2*sind(90+tilt(intx.intAdjacencyMatrix))],...
        'b-', 'LineWidth', 2)
end

k_samp = 2;
trajx = group_traj_x{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
trajy = group_traj_y{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
tilt = group_tilt_dir{4}(1:falloff_traj_ind{4}(fsi(k_samp), k_sub), fsi((k_samp)), k_sub);
plot(trajx, trajy, 'm-', 'LineWidth', 1)
for i_loc = 1:length(k_linsp_1)
    x1 = track_out(1, k_linsp_1(i_loc));
    x2 = track_out(3, k_linsp_2(i_loc));
    y1 = track_out(2, k_linsp_1(i_loc));
    y2 = track_out(4, k_linsp_2(i_loc));
    XY1 = [x1 y1 x2 y2];
    XY2 = [trajx(1:(end-1)), trajy(1:(end-1)), trajx(2:end), trajy(2:end)];
    intx = lineSegmentIntersect(XY1, XY2);
    plot(intx.intMatrixX(intx.intAdjacencyMatrix) + [0 2*cosd(90+tilt(intx.intAdjacencyMatrix))],...
        intx.intMatrixY(intx.intAdjacencyMatrix) + [0 2*sind(90+tilt(intx.intAdjacencyMatrix))],...
        'b-', 'LineWidth', 2)
end