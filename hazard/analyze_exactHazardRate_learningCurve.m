% analyze performance through the distance fall measure (using the exact
% calculation method). Look at Pre-probe and probe epoch performances
% across groups on Probe days and see whether (a) probes are worse than
% pre-probes on days >=3, (b) performance during pre-probe and probe
% increase commensurately wrt one another across days

% David Huberdeau, 06/26/16

%% define globals:
global alignment_error_tolerance
t = .02:.02:1;
alignment_error_tolerance = 0.05;
censor_TH = .95; 
%%

load('exact_track_dist_full_v3.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat') %changed to above,
% 9/27/2017


% load('group_advanced_fall_data.mat')
load('group_track_num.mat')

track_struct_loc = 'TracksStructure_v2.mat';
load(track_struct_loc);

subjects_to_exclude = [8];
subjects_to_exclude_blocks = [4];

%%
figure
subplot(1,2,1); hold on;
trk1 = plotPath(Tracks.fast_challenge_1.Paths);
trk1_l_lens = cumsum(sqrt(diff(trk1(1,:)).^2 + diff(trk1(2,:)).^2));
trk1_r_lens = cumsum(sqrt(diff(trk1(3,:)).^2 + diff(trk1(4,:)).^2));
trk1_l_frac = trk1_l_lens/trk1_l_lens(end);
trk1_r_frac = trk1_r_lens/trk1_r_lens(end);
markers = .1:.1:.9;
for i_mark = 1:length(markers)
    [junkmin, kmin] = min(abs(markers(i_mark) - trk1_l_frac));
    plot(trk1(1, kmin), trk1(2, kmin), 'k.', 'MarkerSize', 15);
    text(trk1(1, kmin), trk1(2, kmin), num2str(markers(i_mark)));
    
    [junkmin, kmin] = min(abs(markers(i_mark) - trk1_r_frac));
    plot(trk1(3, kmin), trk1(4, kmin), 'k.', 'MarkerSize', 15);
    text(trk1(3, kmin), trk1(4, kmin), num2str(markers(i_mark)));
end

subplot(1,2,2); hold on;
trkm = plotPath(Tracks.fast_challenge_1_mirror.Paths);
trkm_l_lens = cumsum(sqrt(diff(trkm(1,:)).^2 + diff(trkm(2,:)).^2));
trkm_r_lens = cumsum(sqrt(diff(trkm(3,:)).^2 + diff(trkm(4,:)).^2));
trkm_l_frac = trkm_l_lens/trkm_l_lens(end);
trkm_r_frac = trkm_r_lens/trkm_r_lens(end);
for i_mark = 1:length(markers)
    [junkmin, kmin] = min(abs(markers(i_mark) - trkm_l_frac));
    plot(trkm(1, kmin), trkm(2, kmin), 'k.', 'MarkerSize', 15);
    text(trkm(1, kmin), trkm(2, kmin), num2str(markers(i_mark)));
    
    [junkmin, kmin] = min(abs(markers(i_mark) - trkm_r_frac));
    plot(trkm(3, kmin), trkm(4, kmin), 'k.', 'MarkerSize', 15);
    text(trkm(3, kmin), trkm(4, kmin), num2str(markers(i_mark)));
end

%% plot overall hazard for training path:

fall_temp_1 = [exact_track_dist_fall{1}; nan(1800, 21)];
fall_temp_1(76:125, :) = nan;

fall_temp_3 = [exact_track_dist_fall{2}; nan(1000, 20)];
fall_temp_3(476:525, :) = nan;

fall_temp_5 = [exact_track_dist_fall{3}; nan(1000, 20)];
fall_temp_5(876:925, :) = nan;

fall_temp_10 = exact_track_dist_fall{4};
fall_temp_10(1876:1925, :) = nan;

pooled_fall_dist = [fall_temp_1,...
    fall_temp_3,...
    fall_temp_5,...
    fall_temp_10];

pooled_haz_func = nan(length(t), size(pooled_fall_dist,2));
for i_sub = 1:size(pooled_haz_func,2)
    try
        l_pooled = hazardRate_nonExtrap(t, pooled_fall_dist(:, i_sub), pooled_fall_dist(:, i_sub) > censor_TH);
        pooled_haz_func(:, i_sub) = l_pooled;
    catch
        pooled_haz_func(:, i_sub) = nan;
    end
end

%%
figure;

% errorfield(t(1:47), sgolayfilt(nanmean(pooled_haz_func(1:47,:), 2), 3, 5), sqrt(nanvar(pooled_haz_func(1:47,:), 0, 2)./sum(~isnan(pooled_haz_func(1,:)))));
errorfield(t(1:47), nanmean(pooled_haz_func(1:47,:), 2), sqrt(nanvar(pooled_haz_func(1:47,:), 0 ,2)./sum(~isnan(pooled_haz_func(1,:)))));

%% plot learning over days for each group:

wind_size = 25;
disc_size = mean(diff(t));

hazardRate_curves = {...
    nan(50, size(exact_track_dist_fall{1}, 1)/wind_size, size(exact_track_dist_fall{1},2)),...
    nan(50, size(exact_track_dist_fall{2}, 1)/wind_size, size(exact_track_dist_fall{2},2)),...
    nan(50, size(exact_track_dist_fall{3}, 1)/wind_size, size(exact_track_dist_fall{3},2)),...
    nan(50, size(exact_track_dist_fall{4}, 1)/wind_size, size(exact_track_dist_fall{4},2))...
    };

mean_hazard_summary = {...
    nan(size(exact_track_dist_fall{1}, 1)/wind_size, size(exact_track_dist_fall{1},2)),...
    nan(size(exact_track_dist_fall{2}, 1)/wind_size, size(exact_track_dist_fall{2},2)),...
    nan(size(exact_track_dist_fall{3}, 1)/wind_size, size(exact_track_dist_fall{3},2)),...
    nan(size(exact_track_dist_fall{4}, 1)/wind_size, size(exact_track_dist_fall{4},2)),...
    };

% range_inds = 13:33; % roi range (indicies: discritized by .02 track fraction)
range_inds = 15:29; %roi range from 0.3 to 0.58
% range_inds = 15:40; %roi range from 0.3 to 0.8
% range_inds = 23:33; %roi range (indicies: discritized by .02 track fraction)
% range_inds = 13:23; %roi range (indicies: discritized by .02 track fraction)

trapezoid_meth = @(dat_mat, range_inds) (disc_size/2)*(dat_mat(range_inds(1),:) +... 
                2*sum(dat_mat(range_inds(2:(end-1)),:)) +... 
                dat_mat(range_inds(end),:));

for i_grp = 1:4
    winds = reshape(1:(wind_size*size(hazardRate_curves{i_grp},2)), wind_size, size(hazardRate_curves{i_grp},2));
    for i_sub = 1:size(exact_track_dist_fall{i_grp},2)
        for i_wind = 1:size(winds,2)
            if sum(isnan(exact_track_dist_fall{i_grp}(winds(:, i_wind), i_sub))) > 5
                % do nothing
            else
%                 l = hazardRate_nonExtrap(t, exact_track_dist_fall{i_grp}(winds(:, i_wind), i_sub),...
%                         exact_track_dist_fall{i_grp}(winds(:, i_wind), i_sub)>censor_TH);

                l = hazardRate_nonExtrap(t, exact_track_dist_fall{i_grp}(winds(:, i_wind), i_sub));

                hazardRate_curves{i_grp}(:, i_wind, i_sub) = l;
            end
%             mean_hazard_summary{i_grp}(i_wind, i_sub) = disc_size*sum(hazardRate_curves{i_grp}(range_inds, i_wind, i_sub));
%             mean_hazard_summary{i_grp}(i_wind, i_sub) = ...
%                 (disc_size/2)*(hazardRate_curves{i_grp}(range_inds(1), i_wind, i_sub) +... 
%                 2*sum(hazardRate_curves{i_grp}(range_inds(2:(end-1)), i_wind, i_sub)) +... 
%                 hazardRate_curves{i_grp}(range_inds(end), i_wind, i_sub));

             mean_hazard_summary{i_grp}(i_wind, i_sub) = trapezoid_meth(hazardRate_curves{i_grp}(:, i_wind, i_sub), range_inds);
        end
    end
end

% exclude subjects:
for i_sub = 1:length(subjects_to_exclude)
    n_windows = size(hazardRate_curves{subjects_to_exclude_blocks(i_sub)}(:,:,subjects_to_exclude(i_sub)), 2);
    hazardRate_curves{subjects_to_exclude_blocks(i_sub)}(:,:,subjects_to_exclude(i_sub)) = nan(50, n_windows);
    mean_hazard_summary{subjects_to_exclude_blocks(i_sub)}(:, subjects_to_exclude(i_sub)) = nan(n_windows,1);
end

figure; hold on;
grp_clr = {'m', 'r', 'b', 'g'};
max_haz = 0;
for i_grp = [4 3 2 1]
    errorfield(1:size(mean_hazard_summary{i_grp},1), ...
        nanmean(mean_hazard_summary{i_grp}, 2),...
        sqrt(nanvar(mean_hazard_summary{i_grp}, 0, 2)./sum(~isnan(mean_hazard_summary{i_grp}),2)),...
        [grp_clr{i_grp}, '-']);
    plot(1:size(mean_hazard_summary{i_grp},1), ...
        nanmean(mean_hazard_summary{i_grp}, 2),...
        [grp_clr{i_grp}, 'o-']);
    
    temp_haz_sd = max(nanmean(mean_hazard_summary{i_grp}, 2) + sqrt(nanvar(mean_hazard_summary{i_grp}, 0, 2)./sum(~isnan(mean_hazard_summary{i_grp}),2)));
    
    if temp_haz_sd > max_haz
        max_haz = temp_haz_sd;
    end
end
for i_day = 1:9
    plot((i_day*[8, 8])+.1, [0, max_haz*1.1], 'Color', [.5 .5 .5]);
end
axis([0 81 0 max_haz*1.1]);
title(['Track range: ', num2str((range_inds(1)/50)*100), '-', num2str((range_inds(end)/50)*100)]);

%% 
hazard_summary_pooled = nan(80, 81);
cancel_trials = [76:125; 476:525; 876:925; 1876:1925];
cancel_trials = (cancel_trials(:, 1:25:50)-1)/25 + 1; % (blocks to cancel)
i_sub = 1;
for i_grp = 1:4
    temp_haz = mean_hazard_summary{i_grp};
    temp_haz(cancel_trials(i_grp, :), :) = nan;
    
    hazard_summary_pooled(1:size(mean_hazard_summary{i_grp},1), ...
        i_sub + (1:size(mean_hazard_summary{i_grp},2)) - 1) = ...
        temp_haz;
    i_sub = i_sub + size(mean_hazard_summary{i_grp},2);
end

figure; hold on;
errorfield(1:75, nanmean(hazard_summary_pooled(1:75,:), 2), sqrt(nanvar(hazard_summary_pooled(1:75, :), 0, 2)./sum(isnan(hazard_summary_pooled(1:75,:)),2)), 'b');
errorfield(78:80, nanmean(hazard_summary_pooled(78:80,:), 2), sqrt(nanvar(hazard_summary_pooled(78:80, :), 0, 2)./sum(isnan(hazard_summary_pooled(78:80,:)),2)), 'b');
for i_day = 1:10
    plot(i_day*[8, 8] + .25, [0 3], 'k-');
end

% export day and across -day window learning data (init)
early_late_inds = reshape([1:8:80; 2:8:80], 20, 1);
early_late_haz = hazard_summary_pooled(early_late_inds, :);
hazard_response = reshape(early_late_haz, 81*20, 1);
csvwrite('hazard_response_init', hazard_response);

% export day and within-day window learning data (early):
early_late_inds = reshape([2:8:80; 8:8:80], 20, 1);
early_late_haz = hazard_summary_pooled(early_late_inds, :);
hazard_response = reshape(early_late_haz, 81*20, 1);
subject_cat = reshape(repmat(1:81, 20, 1), 81*20, 1);
window_cat = reshape(repmat([zeros(1, 81); ones(1,81)], 10, 1), 81*20, 1);
day_cat = repmat(reshape([1:10; 1:10], 20, 1), 81, 1);
hazard_category = [subject_cat, window_cat, day_cat];
csvwrite('hazard_response', hazard_response);
csvwrite('hazard_category', hazard_category);

% export night and across-day window learning data (init):
late_early_inds = reshape([8:8:74; 9:8:74], 18, 1);
late_early_haz = hazard_summary_pooled(late_early_inds, :);
hazard_response_night = reshape(late_early_haz, 81*18, 1);
csvwrite('hazard_response_night_init', hazard_response_night);

% export night and across-day window learning data (early):
late_early_inds = reshape([8:8:74; 10:8:74], 18, 1);
late_early_haz = hazard_summary_pooled(late_early_inds, :);
hazard_response_night = reshape(late_early_haz, 81*18, 1);
subject_cat_night = reshape(repmat(1:81, 18, 1), 81*18, 1);
window_cat_night = reshape(repmat([zeros(1, 81); ones(1,81)], 9, 1), 81*18, 1);
day_cat_night = repmat(reshape([1:9; 1:9], 18, 1), 81, 1);
hazard_category_night = [subject_cat_night, window_cat_night, day_cat_night];
csvwrite('hazard_response_night_early', hazard_response_night);
csvwrite('hazard_category_night', hazard_category_night);

%%

figure;
for i_grp = 1:4
    subplot(2,2,i_grp);
    plot3(repmat((1:50)', 1, size(hazardRate_curves{i_grp},2)), repmat(1:size(hazardRate_curves{i_grp},2), 50, 1), nanmean(hazardRate_curves{i_grp}, 3), grp_clr{i_grp});
end

%% create a plot of hazard within day as a fraction of previous day's asym
n_days = 9;
hazard_summary_day = nan(8, 60*n_days);
pre_day_summary = nan(2, 60*n_days);
i_sub = 1;
day_inds = reshape(9:80, 8, n_days);
prev_asym = 7:8;
for i_day = 1:n_days
    temp = hazard_summary_pooled(day_inds(:, i_day), 22:end)./repmat(nanmean(hazard_summary_pooled(prev_asym, 22:end),1), 8, 1);
    hazard_summary_day(:, i_sub:(i_sub + size(temp,2) - 1)) = temp;
    pre_day_summary(:, i_sub:(i_sub + size(temp,2) - 1)) = hazard_summary_pooled(prev_asym, 22:end)./repmat(nanmean(hazard_summary_pooled(prev_asym, 22:end),1), 2, 1);
    i_sub = i_sub + size(temp,2);
    prev_asym = prev_asym + 8;
end
% hazard_summary_day = reshape(hazard_summary_pooled(9:end, :), 8,  9*81);

figure; hold on
errorbar(1:8, nanmean(hazard_summary_day, 2), sqrt(nanvar(hazard_summary_day, 0, 2)./sum(isnan(hazard_summary_day), 2)), '.-')
plot([0 9], [1 1], 'k-');
ylabel('fraction of previous days asymptotic hazard')
% axis([0, 9, .45, 1.05]);

[h1,p1,c1,d1]=ttest(1 - hazard_summary_day(1,:));
[h2,p2,c2,d2]=ttest(hazard_summary_day(1,:) - hazard_summary_day(8,:));

%% do a cox proportional hazard test for day*probe test

events_all = [exact_track_dist_fall{1}(:); exact_track_dist_fall{2}(:); exact_track_dist_fall{3}(:); exact_track_dist_fall{4}(:)];
c_sub = [reshape(repmat(1:21, 200, 1), 200*21, 1); ...
    reshape(repmat(22:41, 1000, 1), 1000*20, 1); ...
    reshape(repmat(42:61, 1000, 1), 1000*20, 1); ...
    reshape(repmat(62:81, 2000, 1), 2000*20, 1)];
c_grp = [ones(200*21, 1); 2*ones(1000*20, 1); 3*ones(1000*20, 1); 4*ones(2000*20, 1)];
c_day = [ones(200*21, 1);...
    repmat(reshape(repmat(1:5, 200, 1), 200*5, 1), 20, 1);...
    repmat(reshape(repmat(1:5, 200, 1), 200*5, 1), 20, 1);...
    repmat(reshape(repmat(1:10, 200, 1), 200*10, 1), 20, 1)];
c_prb = [reshape([zeros(75, 21); ones(50, 21); zeros(75, 21)], 200*21, 1);...
    reshape([zeros(475, 20); ones(50, 20); zeros(475, 20)], 1000*20, 1);...
    reshape([zeros(875, 20); ones(50, 20); zeros(75, 20)], 1000*20, 1);...
    reshape([zeros(1875, 20); ones(50, 20); zeros(75, 20)], 2000*20, 1)];
events_category = [c_sub, c_grp, c_day, c_prb];
csvwrite('event_times', events_all);
csvwrite('event_category', events_category); %export for coxme in R.
% spoiler: effect of day, effect of probe, significant interaction ** :)

%% do a coxph test for day*window test (within day)
% early_window_inds = [211:260, 411:460, 611:660, 811:860, 1011:1060, 1211:1260, 1411:1460, 1611:1660, 1811:1860];
% asym_window_inds = early_window_inds - 10 + 150; %for window len 50 and excluding first 10 trials
% win_len = 50;%for window len 50 and excluding first 10 trials
early_window_inds = [1:25, 201:225, 401:425, 601:625, 801:825, 1001:1025, 1201:1225, 1401:1425, 1601:1625, 1801:1825];
asym_window_inds = early_window_inds + 175; %for window len 25 and keeping.
win_len = 25;%for window len 25 and keeping first 10 trials
n_days_1 = 1;
n_days_5 = 5;
n_days_10 = 10; % would be 0, 4, and 9 excluding first day
events_within = [reshape(exact_track_dist_fall{2}(early_window_inds(1:(n_days_5*win_len)), :), win_len*20*n_days_5, 1);...
    reshape(exact_track_dist_fall{3}(early_window_inds(1:(n_days_5*win_len)), :), win_len*20*n_days_5, 1); ...
    reshape(exact_track_dist_fall{4}(early_window_inds,:), win_len*20*n_days_10, 1);...
    reshape(exact_track_dist_fall{2}(asym_window_inds(1:(n_days_5*win_len)), :), win_len*20*n_days_5, 1);...
    reshape(exact_track_dist_fall{3}(asym_window_inds(1:(n_days_5*win_len)), :), win_len*20*n_days_5, 1); ...
    reshape(exact_track_dist_fall{4}(asym_window_inds, :), win_len*20*n_days_10, 1)];
c_wind = [zeros(20*n_days_5*win_len*2+20*n_days_10*win_len, 1); ones(20*n_days_5*win_len*2+20*n_days_10*win_len, 1)];
c_sub = [... %no subject 1 bc they only did one day
    reshape(repmat(1:20, win_len*n_days_5, 1), win_len*n_days_5*20, 1); ...
    reshape(repmat(21:40, win_len*n_days_5, 1), win_len*n_days_5*20, 1); ...
    reshape(repmat(41:60, win_len*n_days_10, 1), win_len*n_days_10*20, 1);...
    reshape(repmat(1:20, win_len*n_days_5, 1), win_len*n_days_5*20, 1); ...
    reshape(repmat(21:40, win_len*n_days_5, 1), win_len*n_days_5*20, 1); ...
    reshape(repmat(41:60, win_len*n_days_10, 1), win_len*n_days_10*20, 1)];
c_day = repmat(...
    [repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1); 5*ones(win_len, 1)], 20, 1);...
    repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1); 5*ones(win_len, 1)], 20, 1);...
    repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1); 5*ones(win_len, 1); 6*ones(win_len, 1); 7*ones(win_len, 1); 8*ones(win_len, 1); 9*ones(win_len, 1); 10*ones(win_len, 1)], 20, 1)],...
    2, 1);

events_category_winday = [c_sub, c_day, c_wind];
csvwrite('event_times_winday_25', events_within);
csvwrite('event_category_winday_25', events_category_winday); %export for coxme in R.

%% do a coxph test for day*window test (across days)
% pre_window_inds = [151:200, 351:400, 551:600, 751:800, 951:1000, 1151:1200, 1351:1400, 1551:1600, 1751:1800]; 
% post_window_inds = [211:260, 411:460, 611:660, 811:860, 1011:1060, 1211:1260, 1411:1460, 1611:1660, 1811:1860];
pre_window_inds = [176:200, 376:400, 576:600, 776:800, 976:1000, 1176:1200, 1376:1400, 1576:1600, 1776:1800]; 
post_window_inds = [201:225, 401:425, 601:625, 801:825, 1001:1025, 1201:1225, 1401:1425, 1601:1625, 1801:1825];
win_len = 25;%for window len 25 and keeping first 10 trials
events_within = [reshape(exact_track_dist_fall{2}(pre_window_inds(1:(4*win_len)), :), win_len*20*4, 1);...
    reshape(exact_track_dist_fall{3}(pre_window_inds(1:(4*win_len)), :), win_len*20*4, 1); ...
    reshape(exact_track_dist_fall{4}(pre_window_inds, :), win_len*20*9, 1);...
    reshape(exact_track_dist_fall{2}(post_window_inds(1:(4*win_len)), :), win_len*20*4, 1);...
    reshape(exact_track_dist_fall{3}(post_window_inds(1:(4*win_len)), :), win_len*20*4, 1); ...
    reshape(exact_track_dist_fall{4}(post_window_inds, :), win_len*20*9, 1)];
c_wind = [zeros(20*4*win_len*2+20*9*win_len, 1); ones(20*4*win_len*2+20*9*win_len, 1)];
c_sub = [... %no subject 1 bc they only did one day
    reshape(repmat(1:20, win_len*4, 1), win_len*4*20, 1); ...
    reshape(repmat(21:40, win_len*4, 1), win_len*4*20, 1); ...
    reshape(repmat(41:60, win_len*9, 1), win_len*9*20, 1);...
    reshape(repmat(1:20, win_len*4, 1), win_len*4*20, 1); ...
    reshape(repmat(21:40, win_len*4, 1), win_len*4*20, 1); ...
    reshape(repmat(41:60, win_len*9, 1), win_len*9*20, 1)];
c_day = repmat(...
    [repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1)], 20, 1);...
    repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1)], 20, 1);...
    repmat([ones(win_len, 1); 2*ones(win_len, 1); 3*ones(win_len, 1); 4*ones(win_len, 1); 5*ones(win_len, 1); 6*ones(win_len, 1); 7*ones(win_len, 1); 8*ones(win_len, 1); 9*ones(win_len, 1)], 20, 1)],...
    2, 1);

events_category_winday = [c_sub, c_day, c_wind];
csvwrite('event_times_across', events_within);
csvwrite('event_category_across', events_category_winday); %export for coxme in R.

%% compile and save data for search over blocks for probe/training sim.

events_all = [exact_track_dist_fall{1}(:); exact_track_dist_fall{2}(:); exact_track_dist_fall{3}(:); exact_track_dist_fall{4}(:)];
c_sub = [reshape(repmat(1:21, 200, 1), 200*21, 1); ...
    reshape(repmat(22:41, 1000, 1), 1000*20, 1); ...
    reshape(repmat(42:61, 1000, 1), 1000*20, 1); ...
    reshape(repmat(62:81, 2000, 1), 2000*20, 1)];
c_grp = [ones(200*21, 1); 2*ones(1000*20, 1); 3*ones(1000*20, 1); 4*ones(2000*20, 1)];
c_day = [ones(200*21, 1);...
    repmat(reshape(repmat(1:5, 200, 1), 200*5, 1), 20, 1);...
    repmat(reshape(repmat(1:5, 200, 1), 200*5, 1), 20, 1);...
    repmat(reshape(repmat(1:10, 200, 1), 200*10, 1), 20, 1)];
c_prb = [reshape([zeros(75, 21); ones(50, 21); zeros(75, 21)], 200*21, 1);...
    reshape([zeros(475, 20); ones(50, 20); zeros(475, 20)], 1000*20, 1);...
    reshape([zeros(875, 20); ones(50, 20); zeros(75, 20)], 1000*20, 1);...
    reshape([zeros(1875, 20); ones(50, 20); zeros(75, 20)], 2000*20, 1)];

c_prb2 = [reshape([zeros(80, 21); ones(45, 21); zeros(75, 21)], 200*21, 1);...
    reshape([zeros(480, 20); ones(45, 20); zeros(475, 20)], 1000*20, 1);...
    reshape([zeros(880, 20); ones(45, 20); zeros(75, 20)], 1000*20, 1);...
    reshape([zeros(1880, 20); ones(45, 20); zeros(75, 20)], 2000*20, 1)];

c_pre = [reshape([zeros(25, 21); ones(50, 21); zeros(125, 21)], 200*21, 1);...
    reshape([zeros(425, 20); ones(50, 20); zeros(525, 20)], 1000*20, 1);...
    reshape([zeros(825, 20); ones(50, 20); zeros(125, 20)], 1000*20, 1);...
    reshape([zeros(1825, 20); ones(50, 20); zeros(125, 20)], 2000*20, 1)];
% c_win = [repmat(reshape(repmat(1:4, 50, 1), 200, 1), 21, 1); ...
%     repmat(reshape(repmat(1:20, 50, 1), 1000, 1), 20, 1); ...
%     repmat(reshape(repmat(1:20, 50, 1), 1000, 1), 20, 1); ...
%     repmat(reshape(repmat(1:40, 50, 1), 2000, 1), 20, 1)]; % for windows of 50 trials

c_win = [repmat(reshape(repmat(1:8, 25, 1), 200, 1), 21, 1); ...
    repmat(reshape(repmat(1:40, 25, 1), 1000, 1), 20, 1); ...
    repmat(reshape(repmat(1:40, 25, 1), 1000, 1), 20, 1); ...
    repmat(reshape(repmat(1:80, 25, 1), 2000, 1), 20, 1)]; % for windows of 25 trials

events_category = [c_sub, c_grp, c_day, c_prb, c_prb2, c_pre, c_win];
csvwrite('event_times_window_25_pre_prb2', events_all);
csvwrite('event_category_window_25_pre_prb2', events_category); %export for coxme in R.


%% create test data for survival analysis practice:
n_sub = 10;
n_samp = 50; % per subject per factor
mu1 = .1*randn(1, n_sub) + .25;
mu2 = .1*randn(1, n_sub) + .3;
events1 = nan(length(mu1), n_samp);
events2 = nan(length(mu2), n_samp);
for i_sub = 1:length(mu1)
    events1(i_sub, :) = exprnd(mu1(i_sub), 1, n_samp);
    events2(i_sub, :) = exprnd(mu2(i_sub), 1, n_samp);
end
events1 = reshape(events1, n_samp*n_sub, 1);
events2 = reshape(events2, n_samp*n_sub, 1);

events = [events1; events2];
factor = [zeros(size(events1)); ones(size(events2))];
subject = repmat((1:n_sub)',2*n_samp, 1); 
sample = reshape(repmat(1:50, n_sub, 1), n_samp*n_sub, 1);
category = [subject, factor, [sample; sample]];

% csvwrite('test_survival_events', events);
% csvwrite('test_survival_factor', category);

