
score_fail_grp = cell(1,4);
score_succ_grp = cell(1,4);
slopes_pre_grp = cell(1, 4);
slopes_prb_grp = cell(1, 4);
lm_fail_acrossDays = cell(1, 4);
lm_succ_acrossDays = cell(1, 4);
lm_fail_probeDay = cell(1,4);
lm_succ_probeDay = cell(1,4);
lm_all_probeDay = cell(1,4);
score_fail_by_trial_grp = cell(2, 4); %row1: pre, row2: prb
score_succ_by_trial_grp = cell(2, 4);

for i_grp = 2:4
    %%
    pr_loc_gvn_succ_vPROBE_failures_v5
    score_fail_grp{i_grp} = score;
    lm_fail_acrossDays{i_grp} = lm;
    score_fail_by_trial_grp{1, i_grp} = score_by_trial_pre;
    score_fail_by_trial_grp{2, i_grp} = score_by_trial_prb;
    
    pr_loc_gvn_succ_vPROBE_successes_v5
    score_succ_grp{i_grp} = score;
    lm_succ_acrossDays{i_grp} = lm;
    score_succ_by_trial_grp{1, i_grp} = score_by_trial_pre;
    score_succ_by_trial_grp{2, i_grp} = score_by_trial_prb;
    
    %%
    T_dev = table;
    T_dev.dev = double([T_dev_succ.dev; T_dev_fail.dev]);
    T_dev.len = double([T_dev_succ.len; T_dev_fail.len]);
    T_dev.sub = categorical([T_dev_succ.sub; T_dev_fail.sub]);
    T_dev.succ = [T_dev_succ.succ; categorical(zeros(size(T_dev_fail.succ)))];
    T_dev.prb = categorical([T_dev_succ.prb; T_dev_fail.prb]);
    T_dev.day = categorical([T_dev_succ.day; T_dev_fail.day]);

    prb_day_grp = {1 3 5 10};
%     i_grp = 2;

    lm_all = fitlme(T_dev, 'dev ~ len + succ + len:succ + (1|sub)');
    lm_succ = fitlme(T_dev, 'dev ~ len + prb + prb:len + (1|sub)', 'exclude', T_dev.succ == categorical(0) & T_dev.day ~= categorical(prb_day_grp{i_grp}));
    lm_fail = fitlme(T_dev, 'dev ~ len + prb + prb:len + (1|sub)', 'exclude', T_dev.succ == categorical(1) & T_dev.day ~= categorical(prb_day_grp{i_grp}));

    lm_fail_probeDay{i_grp} = lm_fail;
    lm_succ_probeDay{i_grp} = lm_succ;
    lm_all_probeDay{i_grp} = lm_all;
    
    slopes_pre = nan(2, n_subs);
    slopes_prb = nan(2, n_subs);

    for i_sub = setdiff(1:n_subs, REMOVE_SUB)

        sub_temp = T_dev.sub;
        len_temp = T_dev.len;
        dev_temp = T_dev.dev;
        suc_temp = T_dev.succ;
        prb_temp = T_dev.prb;
        day_temp = T_dev.day;

        T_temp = table;
        T_temp.len = len_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));
        T_temp.sub = sub_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));
        T_temp.dev = dev_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));

        try
            lm_temp = fitlme(T_temp, 'dev ~ len');
            slopes_pre(1, i_sub) = lm_temp.Coefficients.Estimate(2);
        catch
            slopes_pre(1, i_sub) = nan;
        end

        T_temp = table;
        T_temp.len = len_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));
        T_temp.sub = sub_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));
        T_temp.dev = dev_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(0));

        try
            lm_temp = fitlme(T_temp, 'dev ~ len');
            slopes_pre(2, i_sub) = lm_temp.Coefficients.Estimate(2);
        catch
            slopes_pre(2, i_sub) = nan;
        end

        T_temp = table;
        T_temp.len = len_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));
        T_temp.sub = sub_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));
        T_temp.dev = dev_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(1) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));

        try
            lm_temp = fitlme(T_temp, 'dev ~ len');
            slopes_prb(1, i_sub) = lm_temp.Coefficients.Estimate(2);
        catch
            slopes_prb(1, i_sub) = nan;
        end

        T_temp = table;
        T_temp.len = len_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));
        T_temp.sub = sub_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));
        T_temp.dev = dev_temp(sub_temp == categorical(i_sub) & suc_temp == categorical(0) & day_temp == categorical(prb_day_grp{i_grp}) & prb_temp == categorical(1));

        try
            lm_temp = fitlme(T_temp, 'dev ~ len');
            slopes_prb(2, i_sub) = lm_temp.Coefficients.Estimate(2);
        catch
            slopes_prb(2, i_sub) = nan;
        end
    end
    slopes_pre_grp{i_grp} = slopes_pre;
    slopes_prb_grp{i_grp} = slopes_prb;
    close all;
    
    
end
%%
% slope analysis

slopes_pre = [slopes_pre_grp{2};slopes_pre_grp{3};slopes_pre_grp{4}];
slopes_prb = [slopes_prb_grp{2};slopes_prb_grp{3};slopes_prb_grp{4}];

figure; hold on;
errorbar([1 2], nanmean([slopes_pre(2,:)', slopes_prb(2,:)']), nanstd([slopes_pre(2,:)', slopes_prb(2,:)'])./sqrt(17), 'r.-')
errorbar([1 2], nanmean([slopes_pre(1,:)', slopes_prb(1,:)']), nanstd([slopes_pre(1,:)', slopes_prb(1,:)'])./sqrt(17), 'b.-')
title('Slope')

figure;
for i_grp = 2:4
    slopes_pre = [slopes_pre_grp{i_grp}];
    slopes_prb = [slopes_prb_grp{i_grp}];

    subplot(1,4, i_grp); hold on;
    errorbar([1 2], nanmean([slopes_pre(2,:)', slopes_prb(2,:)']), nanstd([slopes_pre(2,:)', slopes_prb(2,:)'])./sqrt(17), 'r.-')
    errorbar([1 2], nanmean([slopes_pre(1,:)', slopes_prb(1,:)']), nanstd([slopes_pre(1,:)', slopes_prb(1,:)'])./sqrt(17), 'b.-')
    title('Slope')
    axis([0 3 -4 16])
end

lm_slope = cell(1,4);
for i_grp = 2:4
    T_slope = table;
    T_slope.slope = [slopes_pre_grp{i_grp}(1,:)'; slopes_pre_grp{i_grp}(2,:)';...
        slopes_prb_grp{i_grp}(1,:)'; slopes_prb_grp{i_grp}(2,:)'];
    T_slope.subject = repmat((1:20)', 4, 1);
    T_slope.probe = [zeros(40, 1); ones(40, 1)];
    T_slope.success = [ones(20,1); zeros(20,1); ones(20,1); zeros(20,1)];
    
    lm_slope{i_grp} = fitlme(T_slope, 'slope ~ probe + success + probe:success + (1|subject)');
end

%%
% mean deviation analysis

score_succ_all = [[score_succ_grp{2}(:, 1:5), nan(20, 1),score_succ_grp{2}(:, 7:10), nan(20, 10)];...
    [score_succ_grp{3}(:, setdiff(1:10, 10)), nan(20, 11)];...
    score_succ_grp{4}(:, setdiff(1:20, 20)), nan(20, 1)];
score_succ_prbDay = [score_succ_grp{2}(:, 5:6); score_succ_grp{3}(:, 9:10); score_succ_grp{4}(:, 19:20)];


score_fail_all = [[score_fail_grp{2}(:, 1:5), nan(20, 1),score_fail_grp{2}(:, 7:10), nan(20, 10)];...
    [score_fail_grp{3}(:, setdiff(1:10, 10)), nan(20, 11)];...
    score_fail_grp{4}(:, setdiff(1:20, 20)), nan(20, 1)];
score_fail_prbDay = [score_fail_grp{2}(:, 5:6); score_fail_grp{3}(:, 9:10); score_fail_grp{4}(:, 19:20)];


[h_score_s, p_score_s, b_score_s, c_score_s] = ttest(score_succ_prbDay(:, 1) - score_succ_prbDay(:, 2));
[h_score_f, p_score_f, b_score_f, c_score_f] = ttest(score_fail_prbDay(:, 1) - score_fail_prbDay(:, 2));

%%
figure; hold on;
plot([1 2], [score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)], '-', 'Color', [.5 .5 .5]);
errorbar([1 2], nanmean([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)]),...
    nanstd([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)])./sqrt(sum(~isnan([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)]))),...
    'b.-');
axis([.8 2.2 .2 2])
title('Successes');
figure; hold on;
plot([1 2], [score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)], '-', 'Color', [.5 .5 .5]);
errorbar([1 2], nanmean([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)]),...
    nanstd([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)])./sqrt(sum(~isnan([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)]))),...
    'r.-');
axis([.8 2.2 .2 2])
title('Failures');

%%
dt_pre_all = nan(60, 10);
dt_prb_all = nan(60, 10);
k_subj = 1;
for i_grp = 2:4
    dt_pre = nan(20, grp_days(i_grp));
    dt_prb = nan(20, grp_days(i_grp));
    for i_day = 1:grp_days(i_grp)
        for i_sub = 1:20
            valid_inds = exact_track_dist_fall{i_grp}(pre_wind{i_day}, i_sub) > .5;
            dt_pre(i_sub, i_day) = nanmean(exact_track_dist_fall{i_grp}(pre_wind{i_day}(valid_inds), i_sub), 1)';
            valid_inds = exact_track_dist_fall{i_grp}(probe_wind{i_day}, i_sub) > .5;
            dt_prb(i_sub, i_day) = nanmean(exact_track_dist_fall{i_grp}(probe_wind{i_day}(valid_inds), i_sub), 1)';
        end
    end
    dt_pre_all(k_subj - 1 + (1:20), 1:grp_days(i_grp)) = dt_pre;
    dt_prb_all(k_subj - 1 + (1:20), 1:grp_days(i_grp)) = dt_prb;
    k_subj = k_subj + 20;
end
dt_succ_all = reshape([dt_pre_all; dt_prb_all], 60, 20);


figure; hold on;

succ_x = score_succ_all(:);
succ_y = dt_succ_all(:);

throw_inds = succ_x > 1.8;
succ_x(throw_inds) = nan;
succ_y(throw_inds) = nan;
lin_succ = fitlm(succ_x, succ_y);

plot(score_succ_all(:), dt_succ_all(:), 'k.');
plot([min(succ_x), max(succ_x)], lin_succ.Coefficients.Estimate(1) + lin_succ.Coefficients.Estimate(2)*[min(succ_x), max(succ_x)], 'b-');

dt_succ_prbDay = [dt_succ_all(1:20, 5:6); dt_succ_all(21:40, 9:10); dt_succ_all(41:60, 19:20)];
pre_resid = nan(size(score_succ_prbDay,1), 1);
prb_resid = nan(size(score_succ_prbDay,1), 1);
for i_resid = 1:length(pre_resid)
    pre_resid(i_resid) = score_succ_prbDay(i_resid,1) - (lin_succ.Coefficients.Estimate(1) + lin_succ.Coefficients.Estimate(2)*dt_succ_prbDay(i_resid,1));
    prb_resid(i_resid) = score_succ_prbDay(i_resid,2) - (lin_succ.Coefficients.Estimate(1) + lin_succ.Coefficients.Estimate(2)*dt_succ_prbDay(i_resid,2));
end
plot(score_succ_prbDay(:,1), dt_succ_prbDay(:,1), 'g.');
plot(score_succ_prbDay(:,2), dt_succ_prbDay(:,2), 'r.');

figure;
hold on;
plot([1 2], [pre_resid, prb_resid], '-', 'Color', [.5 .5 .5]);
errorbar([1, 2], nanmean([pre_resid, prb_resid]), nanstd([pre_resid, prb_resid])./sqrt(sum(~isnan([pre_resid, prb_resid]))), 'k-', 'LineWidth', 3);
[h_resid, p_resid, a_resid, b_resid] = ttest(prb_resid - pre_resid);


figure; hold on;
% plot([1 2], [score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)], '-', 'Color', [.5 .5 .5]);
errorbar([1 2], nanmean([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)]),...
    nanstd([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)])./sqrt(sum(~isnan([score_succ_prbDay(:, 1), score_succ_prbDay(:, 2)]))),...
    'b.-');
% title('Successes');
% figure; hold on;
% plot([1 2], [score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)], '-', 'Color', [.5 .5 .5]);
errorbar([1 2], nanmean([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)]),...
    nanstd([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)])./sqrt(sum(~isnan([score_fail_prbDay(:, 1), score_fail_prbDay(:, 2)]))),...
    'r.-');
title('Deviation interaction');


lm_dev = cell(1,4);
sub_grp_k = {[], 1:20, 21:40, 41:60};
for i_grp = 2:4
    T_mdev = table;
    T_mdev.dev = [score_succ_prbDay(sub_grp_k{i_grp}, 1); score_fail_prbDay(sub_grp_k{i_grp}, 1);...
        score_succ_prbDay(sub_grp_k{i_grp}, 2); score_fail_prbDay(sub_grp_k{i_grp}, 2)];
    T_mdev.subject = repmat((1:20)', 4, 1);
    T_mdev.probe = [zeros(40, 1); ones(40, 1)];
    T_mdev.success = [ones(20,1); zeros(20,1); ones(20,1); zeros(20,1)];
    
    lm_dev{i_grp} = fitlme(T_mdev, 'dev ~ probe + success + probe:success + (1|subject)');
end

save score_by_trial_grp score_succ_by_trial_grp score_fail_by_trial_grp

%%
pre_wind_wind = [[.375 1.375, 2.375, 3.375 4.375], 5+[.375 1.375, 2.375, 3.375 4.375]];
prb_wind_wind = [[.625 1.625, 2.625, 3.625 4.625], 5+[.625 1.625, 2.625, 3.625 4.625]];
x = reshape([pre_wind_wind(1:grp_days(i_grp)); prb_wind_wind(1:grp_days(i_grp))], 2*grp_days(i_grp), 1);

%learning curve with probes:
figure; %subplot(2,1,1);
errorbar(x, nanmean(score_succ_all,1), sqrt(nanvar(score_succ_all)./sum(~isnan(score_succ_all))), 'b.-')
hold on
errorbar(x(6), nanmean(score_succ_grp{2}(:, 6)), sqrt(nanvar(score_succ_grp{2}(:, 6))./sum(~isnan(score_succ_grp{2}(:, 6)))), 'rs')
errorbar(x(10), nanmean(score_succ_grp{3}(:, 10)), sqrt(nanvar(score_succ_grp{3}(:, 10))./sum(~isnan(score_succ_grp{3}(:, 10)))), 'rs')
errorbar(x(20), nanmean(score_succ_grp{4}(:, 20)), sqrt(nanvar(score_succ_grp{4}(:, 20))./sum(~isnan(score_succ_grp{4}(:, 20)))), 'rs')

% subplot(2,1,2);
errorbar(x, nanmean(score_fail_all,1), sqrt(nanvar(score_fail_all)./sum(~isnan(score_fail_all))), 'g.-')
hold on
errorbar(x(6), nanmean(score_fail_grp{2}(:, 6)), sqrt(nanvar(score_fail_grp{2}(:, 6))./sum(~isnan(score_fail_grp{2}(:, 6)))), 'rs')
errorbar(x(10), nanmean(score_fail_grp{3}(:, 10)), sqrt(nanvar(score_fail_grp{3}(:, 10))./sum(~isnan(score_fail_grp{3}(:, 10)))), 'rs')
errorbar(x(20), nanmean(score_fail_grp{4}(:, 20)), sqrt(nanvar(score_fail_grp{4}(:, 20))./sum(~isnan(score_fail_grp{4}(:, 20)))), 'rs')