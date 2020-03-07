%% plot response to probe, using pre-probe as baseline.
pre_probe_inds = {[26:75], [426:475], [826:875], [1826:1875]};
probe_inds = {[76:125], [476:525], [876:925], [1876:1925]};

pre_wind_scores = {nan(50, 20), nan(50, 20), nan(50, 20), nan(50, 20)};
prb_wind_scores = {nan(50, 20), nan(50, 20), nan(50, 20), nan(50, 20)};

group_colors = {'y', 'r', 'b', 'g'};

for i_grp = 1:4
    score_avg_grp = reshape(nanmean(score_s{i_grp}, 2), ...
        size(score_s{i_grp},1), size(score_s{i_grp},3));
    
    pre_score = score_avg_grp(pre_probe_inds{i_grp}, :);
    prb_score = score_avg_grp(probe_inds{i_grp}, :);
    
    pre_wind_scores{i_grp} = pre_score;
    prb_wind_scores{i_grp} = prb_score;
end


%% plot each day's learning curve, using prior-day's asymptote as baseline
grp_days = [1, 3, 5, 10];
day_baseline_inds = repmat((100:200:1900)', 1, 50) + repmat(51:100, 10, 1);
all_inds = 1:2000;
all_inds_excl = all_inds;
for i_grp = 1:4
    all_inds_excl(probe_inds{i_grp}) = nan;
end

day_nonprobe_inds = {...
    reshape(all_inds_excl(1:200), 200, grp_days(1))' ...
    reshape(all_inds_excl(1:600),200, grp_days(2))', ...
    [reshape(all_inds(1:800), 200, grp_days(3)-1)'; ...
    reshape(all_inds_excl(801:1000), 200, 1)'],...
    [reshape(all_inds(1:1800), 200, grp_days(4)-1)'; ...
    reshape(all_inds_excl(1801:2000), 200, 1)']};

baseline_score_grp = {nan(50, 20), nan(50, 20, grp_days(2)),...
    nan(50, 20, grp_days(3)), nan(50, 20, grp_days(4))};
day_score_grp = {nan(200, 20), nan(200, 20, grp_days(2)),...
    nan(200, 20, grp_days(3)), nan(200, 20, grp_days(4))};

for i_grp = 2:4
    for i_day = 1:grp_days(i_grp)
        score_avg_grp = reshape(nanmean(score_s{i_grp}, 2), ...
            size(score_s{i_grp},1), size(score_s{i_grp},3));
    
        score_baseline_day = score_avg_grp(day_baseline_inds(i_day, :), :);
        baseline_score_grp{i_grp}(:, :, i_day) = score_baseline_day;
        
        all_inds = day_nonprobe_inds{i_grp}(i_day, :);
        valid_inds = all_inds(~isnan(all_inds));
        score_learning_day = score_avg_grp(valid_inds, :);
        day_score_grp{i_grp}(valid_inds - 200*(i_day-1), :, i_day) = score_learning_day;
    end
end
    
grp_inds = {1:21, 21+(1:20), 21+20+(1:20), 21+20+20+(1:20)};
baseline_scores = nan(50, 81);
day_scores = nan(200, 81);
for i_grp = 2:4
    baseline_scores(:, grp_inds{i_grp}) = nanmean(baseline_score_grp{i_grp}, 3);
    day_scores(:, grp_inds{i_grp}) = nanmean(day_score_grp{i_grp}, 3);
end

%%
% figure; hold on;
% plot(1:5:50, nanmean(mafilt(baseline_scores, 5, 1),2), 'Color', [.5 .5 .5]);
all_baseline_scores = [baseline_scores, pre_wind_scores{1}, ...
    pre_wind_scores{2}, pre_wind_scores{3}, pre_wind_scores{4}];

baseline_curve = nanmean(mafilt(all_baseline_scores, 5, 1), 2);

figure; hold on;
plot(-50:5:-1,...
    baseline_curve/mean(baseline_curve),...
    'Color', [.5 .5 .5]);
plot(1:5:200, nanmean(mafilt(day_scores, 5, 1)/nanmean(nanmean(baseline_scores,1),2),2), 'Color', [.5 .5 .5]);
% plot([-50 200], [1 1], 'k-')
for i_grp = 1:4
%     subplot(1,4,i_grp);
%     hold on;
%     plot(1:5:50, nanmean(mafilt(pre_wind_scores{i_grp}, 5, 1),2), 'Color', [.5 .5 .5]);
    plot(1:5:50, nanmean(mafilt(prb_wind_scores{i_grp}, 5, 1)/...
        nanmean(nanmean(pre_wind_scores{i_grp},1),2),2), ...
        group_colors{i_grp});
end