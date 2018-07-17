% analyze performance through the distance fall measure (using the exact
% calculation method). Look at Pre-probe and probe epoch performances
% across groups on Probe days and see whether (a) probes are worse than
% pre-probes on days >=3, (b) performance during pre-probe and probe
% increase commensurately wrt one another across days

% David Huberdeau, 06/26/16

% load('exact_track_dist_fall.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat')
load('exact_track_dist_full_v3');

preProbe_trials = {26:75, 426:475, 826:875, 1826:1875};
probe_trials = {76:125, 476:525, 876:925, 1876:1925};
probe_trials_xbg = {86:125, 486:525, 886:925, 1886:1925};

n_subs_valid = 68;
%%

binsize = 5;
learn_curve_mean = nan(ceil(2000/binsize), 1);
learn_curve_sd = nan(ceil(2000/binsize), 1);
probe_curve_mean = nan(ceil(2000/binsize), 1);
probe_curve_sd = nan(ceil(2000/binsize), 1);

train_only_trials0 = [1:75, 126:200];
train_only_bins0 = train_only_trials0(5:5:length(train_only_trials0))/5;
train_only_trials1 = [201:475, 526:875, 926:1000];
train_only_bins1 = train_only_trials1(5:5:length(train_only_trials1))/5;
train_only_trials2 = [1001:1875, 1926:2000];
train_only_bins2 = train_only_trials2(5:5:length(train_only_trials2))/5;

learn_curve_mean(train_only_bins0) = nanmean(mafilt([exact_track_dist_fall{1}(train_only_trials0,:), ...
    exact_track_dist_fall{2}(train_only_trials0,:),...
    exact_track_dist_fall{3}(train_only_trials0,:)...
    exact_track_dist_fall{4}(train_only_trials0,:)], 5), 2);
learn_curve_sd(train_only_bins0) = sqrt(nanvar(mafilt([exact_track_dist_fall{1}(train_only_trials0,:), ...
    exact_track_dist_fall{2}(train_only_trials0,:),...
    exact_track_dist_fall{3}(train_only_trials0,:)...
    exact_track_dist_fall{4}(train_only_trials0,:)], 5), 0, 2)/n_subs_valid);

learn_curve_mean(train_only_bins1) = nanmean([mafilt(exact_track_dist_fall{2}(train_only_trials1,:), 5),...
    mafilt(exact_track_dist_fall{3}(train_only_trials1,:), 5),...
    mafilt(exact_track_dist_fall{4}(train_only_trials1,:), 5)] ,2);
learn_curve_sd(train_only_bins1) = sqrt(nanvar([mafilt(exact_track_dist_fall{2}(train_only_trials1,:), 5),...
    mafilt(exact_track_dist_fall{3}(train_only_trials1,:), 5),...
    mafilt(exact_track_dist_fall{4}(train_only_trials1,:), 5)], 0, 2)/n_subs_valid);

learn_curve_mean(train_only_bins2) = nanmean(mafilt(exact_track_dist_fall{4}(train_only_trials2,:), 5), 2);
learn_curve_sd(train_only_bins2) = sqrt(nanvar(mafilt(exact_track_dist_fall{4}(train_only_trials2,:), 5), 0, 2)/n_subs_valid);

all_trials = 1:2000;
all_bins = all_trials(5:5:length(all_trials))/5;
grp_probe_trials = {76:126, 476:525, 876:925, 1876:1925};
grp_probe_bins = cell(1,4);
set_compare = {[1 2 3 4], [2 3 4], [2 3 4], []};
for i_grp = [1 2 3 4]
    grp_probe_bins{i_grp} = grp_probe_trials{i_grp}(5:5:length(grp_probe_trials{i_grp}))/5;
    probe_curve_mean(grp_probe_bins{i_grp}) = nanmean(mafilt(exact_track_dist_fall{i_grp}(grp_probe_trials{i_grp},:), 5), 2);
    probe_curve_sd(grp_probe_bins{i_grp}) = sqrt(nanvar(mafilt(exact_track_dist_fall{i_grp}(grp_probe_trials{i_grp},:), 5), 0, 2)/n_subs_valid);
    train_grps = setdiff(set_compare{i_grp}, i_grp);
    temp_train_curve = nan(10, 81); i_sub_temp = 1;
    if ~isempty(train_grps)
        for i_train_grp = train_grps
            temp_train_curve(:, i_sub_temp:(i_sub_temp + size(exact_track_dist_fall{i_train_grp},2) - 1)) = mafilt(exact_track_dist_fall{i_train_grp}(grp_probe_trials{i_grp},:),5);
            i_sub_temp = i_sub_temp + size(exact_track_dist_fall{i_train_grp},2);
        end
    end
    learn_curve_mean(grp_probe_bins{i_grp}) = nanmean(temp_train_curve, 2);
    learn_curve_sd(grp_probe_bins{i_grp}) = sqrt(nanvar(temp_train_curve, 0 ,2)/n_subs_valid);
end

figure;
errorfield(1:5:2000, learn_curve_mean, learn_curve_sd, 'k'); hold on
errorfield(1:5:2000, probe_curve_mean, probe_curve_sd, 'r'); hold on

%%
grp_clr = {'m', 'r', 'b', 'k'};
binsize = 5;
figure; hold on;
for i_grp = [4 3 2 1]
    ytemp = nanmean(mafilt(exact_track_dist_fall{i_grp}, binsize), 2);
    xtemp = binsize*(1:length(ytemp))-binsize/2;
    ztemp = sqrt(nanvar(mafilt(exact_track_dist_fall{i_grp}, binsize), 0, 2)./sum(~isnan(mafilt(exact_track_dist_fall{i_grp}, binsize)),2));
    errorfield(xtemp, ytemp, ztemp, grp_clr{i_grp});
end
for i_day = 1:10
    plot(i_day*[200.5 200.5], [0 1], 'k');
end
plot([75.5 75.5], [0 1], 'k--'); plot([125.5 125.5], [0 1], 'k--');
plot([475.5 475.5], [0 1], 'k--'); plot([525.5 525.5], [0 1], 'k--');
plot([875.5 875.5], [0 1], 'k--'); plot([925.5 925.5], [0 1], 'k--');
plot([1875.5 1875.5], [0 1], 'k--'); plot([1925.5 1925.5], [0 1], 'k--');
axis([0 2000 .2 .8])

%%
avgDist_pre = nan(size(exact_track_dist_fall{1},2), length(exact_track_dist_fall));
avgDist_prb = nan(size(exact_track_dist_fall{1},2), length(exact_track_dist_fall));
avgDist_prb_xbg = nan(size(exact_track_dist_fall{1},2), length(exact_track_dist_fall));
avgDist_prb_bg = nan(size(exact_track_dist_fall{1},2), length(exact_track_dist_fall));

for i_grp = 1:length(exact_track_dist_fall)
    avgDist_pre(1:size(exact_track_dist_fall{i_grp},2), i_grp) = ...
        nanmean(exact_track_dist_fall{i_grp}(preProbe_trials{i_grp}, :), 1)';
    
    avgDist_prb(1:size(exact_track_dist_fall{i_grp},2), i_grp) = ...
        nanmean(exact_track_dist_fall{i_grp}(probe_trials{i_grp}, :), 1)';
    
    avgDist_prb_xbg(1:size(exact_track_dist_fall{i_grp},2), i_grp) = ...
        nanmean(exact_track_dist_fall{i_grp}(probe_trials_xbg{i_grp}, :), 1)';
    
    avgDist_prb_bg(1:size(exact_track_dist_fall{i_grp},2), i_grp) = ...
        nanmean(exact_track_dist_fall{i_grp}(setdiff(probe_trials{i_grp}, probe_trials_xbg{i_grp}), :), 1)';
end

%% view group data on probe day:
figure; %hold on
subplot(2,2,1); hold on;
temp=nanmean(mafilt(exact_track_dist_fall{1}(1:200,:),5), 2);
temp2=sqrt(nanvar(mafilt(exact_track_dist_fall{1}(1:200,:),5),0, 2)./sum(~isnan(mafilt(exact_track_dist_fall{1}(1:200,:),5)),2));
errorfield([2.5:5:200], temp, temp2, 'm')

temp_data = [mafilt(exact_track_dist_fall{2}(1:200,:),5),...
    mafilt(exact_track_dist_fall{3}(1:200,:),5),...
    mafilt(exact_track_dist_fall{4}(1:200,:),5)];
temp = nanmean(temp_data, 2);
temp2 = sqrt(nanvar(temp_data, 0, 2)./sum(~isnan(temp_data),2));
errorfield(2.5:5:200, temp, temp2, 'k')
plot([75.5 75.5], mean(temp)+[-0.25 0.25], 'k--')
plot([125.5 125.5], mean(temp)+[-0.25 0.25], 'k--')
% axis([0 200 mean(temp)+[-0.25 0.25]])
axis([0 200 0.2 0.9])

subplot(2,2,2); hold on;
temp=nanmean(mafilt(exact_track_dist_fall{2}(401:600,:),5), 2);
temp2=sqrt(nanvar(mafilt(exact_track_dist_fall{2}(401:600,:),5),0, 2)./sum(~isnan(mafilt(exact_track_dist_fall{2}(401:600,:),5)),2));
errorfield([2.5:5:200], temp, temp2, 'r')

temp_data = [mafilt(exact_track_dist_fall{3}(401:600,:),5),...
    mafilt(exact_track_dist_fall{4}(401:600,:),5)];
temp = nanmean(temp_data, 2);
temp2 = sqrt(nanvar(temp_data, 0, 2)./sum(~isnan(temp_data),2));
errorfield(2.5:5:200, temp, temp2, 'k')
plot([75.5 75.5], mean(temp)+[-0.25 0.25], 'k--')
plot([125.5 125.5], mean(temp)+[-0.25 0.25], 'k--')
% axis([0 200 mean(temp)+[-0.25 0.25]])
axis([0 200 0.2 0.9])

subplot(2,2,3); hold on;
temp=[...nanmean(exact_track_dist_fall{3}(801:810,:),2); ...
    nanmean(mafilt(exact_track_dist_fall{3}(801:1000,:),5), 2)];
temp2=sqrt([...nanvar(exact_track_dist_fall{3}(801:810,:),0,2)./sum(~isnan(exact_track_dist_fall{3}(801:810,:)),2); ...
    nanvar(mafilt(exact_track_dist_fall{3}(801:1000,:),5),0, 2)./sum(~isnan(mafilt(exact_track_dist_fall{3}(801:1000,:),5)),2)]);
errorfield([2.5:5:200], temp, temp2, 'b');

temp_data = [mafilt(exact_track_dist_fall{2}(801:1000,:),5),...
    mafilt(exact_track_dist_fall{4}(801:1000,:),5)];
temp = nanmean(temp_data, 2);
temp2 = sqrt(nanvar(temp_data, 0, 2)./sum(~isnan(temp_data),2));
errorfield(2.5:5:200, temp, temp2, 'k')
plot([75.5 75.5], mean(temp)+[-0.25 0.25], 'k--')
plot([125.5 125.5], mean(temp)+[-0.25 0.25], 'k--')
% axis([0 200 mean(temp)+[-0.25 0.25]])
axis([0 200 0.2 0.9])

subplot(2,2,4); hold on;
temp = nanmean(mafilt(exact_track_dist_fall{4}(1801:2000,:),5), 2);
temp2=sqrt(nanvar(mafilt(exact_track_dist_fall{4}(1801:2000,:),5),0, 2)./sum(~isnan(mafilt(exact_track_dist_fall{4}(1801:2000,:),5)),2));
errorfield([2.5:5:200], temp, temp2, 'g')

plot([75.5 75.5], mean(temp)+[-0.25 0.25], 'k--')
plot([125.5 125.5], mean(temp)+[-0.25 0.25], 'k--')
% axis([0 200 mean(temp)+[-0.25 0.25]])
axis([0 200 0.2 0.9])

% grpmean_avgDist_pre = nanmean(avgDist_pre);
% grpmean_avgDist_prb = nanmean(avgDist_prb);
% grpmean_avgDist_prb_xbg = nanmean(avgDist_prb_xbg);
% grpmean_avgDist_prb_bg = nanmean(avgDist_prb_bg);
% 
% plot([26 75], [grpmean_avgDist_pre(1) grpmean_avgDist_pre(1)], 'm-')
% plot([26 75], [grpmean_avgDist_pre(2) grpmean_avgDist_pre(2)], 'r-')
% plot([26 75], [grpmean_avgDist_pre(3) grpmean_avgDist_pre(3)], 'b-')
% plot([26 75], [grpmean_avgDist_pre(4) grpmean_avgDist_pre(4)], 'k-')
% 
% plot([76 125], [grpmean_avgDist_prb(1) grpmean_avgDist_prb(1)], 'm-')
% plot([76 125], [grpmean_avgDist_prb(2) grpmean_avgDist_prb(2)], 'r-')
% plot([76 125], [grpmean_avgDist_prb(3) grpmean_avgDist_prb(3)], 'b-')
% plot([76 125], [grpmean_avgDist_prb(4) grpmean_avgDist_prb(4)], 'k-')
% 
% plot([81 125], [grpmean_avgDist_prb_xbg(1) grpmean_avgDist_prb_xbg(1)], 'm--')
% plot([81 125], [grpmean_avgDist_prb_xbg(2) grpmean_avgDist_prb_xbg(2)], 'r--')
% plot([81 125], [grpmean_avgDist_prb_xbg(3) grpmean_avgDist_prb_xbg(3)], 'b--')
% plot([81 125], [grpmean_avgDist_prb_xbg(4) grpmean_avgDist_prb_xbg(4)], 'k--')
% 
% plot([76 81], [grpmean_avgDist_prb_bg(1) grpmean_avgDist_prb_bg(1)], 'm.-')
% plot([76 81], [grpmean_avgDist_prb_bg(2) grpmean_avgDist_prb_bg(2)], 'r.-')
% plot([76 81], [grpmean_avgDist_prb_bg(3) grpmean_avgDist_prb_bg(3)], 'b.-')
% plot([76 81], [grpmean_avgDist_prb_bg(4) grpmean_avgDist_prb_bg(4)], 'k.-')

%% view bar plots of data and differences:

figure;
subplot(231); hold on
bar(1, nanmean(avgDist_pre(:,1)), 'm');
bar(2, nanmean(avgDist_pre(:,2)), 'r');
bar(3, nanmean(avgDist_pre(:,3)), 'b');
bar(4, nanmean(avgDist_pre(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_pre), sqrt(nanvar(avgDist_pre, 0, 1)./sum(~isnan(avgDist_pre))), 'k.')
title('Pre-probe: dist');
axis([0 5 0 1])

subplot(232); hold on
bar(1, nanmean(avgDist_prb(:,1)), 'm');
bar(2, nanmean(avgDist_prb(:,2)), 'r');
bar(3, nanmean(avgDist_prb(:,3)), 'b');
bar(4, nanmean(avgDist_prb(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_prb), sqrt(nanvar(avgDist_prb, 0, 1)./sum(~isnan(avgDist_prb))), 'k.')
title('Probe: dist');
axis([0 5 0 1])

subplot(233); hold on
bar(1, nanmean(avgDist_prb_xbg(:,1)), 'm');
bar(2, nanmean(avgDist_prb_xbg(:,2)), 'r');
bar(3, nanmean(avgDist_prb_xbg(:,3)), 'b');
bar(4, nanmean(avgDist_prb_xbg(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_prb_xbg), sqrt(nanvar(avgDist_prb_xbg, 0, 1)./sum(~isnan(avgDist_prb_xbg))), 'k.')
title('Probe Xbegin: dist');
axis([0 5 0 1])


subplot(234); hold on
bar(1, nanmean(avgDist_prb(:,1) - avgDist_pre(:,1)), 'm');
bar(2, nanmean(avgDist_prb(:,2) - avgDist_pre(:,2)), 'r');
bar(3, nanmean(avgDist_prb(:,3) - avgDist_pre(:,3)), 'b');
bar(4, nanmean(avgDist_prb(:,4) - avgDist_pre(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_prb - avgDist_pre), ...
    sqrt(nanvar(avgDist_prb - avgDist_pre, 0, 1)./sum(~isnan(avgDist_prb - avgDist_pre))), 'k.')
title('Probe - Pre: diff');
axis([0 5 -.25 .05])


subplot(235); hold on
bar(1, nanmean(avgDist_prb_bg(:,1) - avgDist_pre(:,1)), 'm');
bar(2, nanmean(avgDist_prb_bg(:,2) - avgDist_pre(:,2)), 'r');
bar(3, nanmean(avgDist_prb_bg(:,3) - avgDist_pre(:,3)), 'b');
bar(4, nanmean(avgDist_prb_bg(:,4) - avgDist_pre(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_prb_bg - avgDist_pre), ...
    sqrt(nanvar(avgDist_prb_bg - avgDist_pre, 0, 1)./sum(~isnan(avgDist_prb_bg - avgDist_pre))), 'k.')
title('Probe begin - Pre: diff');
axis([0 5 -.25 .05])

subplot(236); hold on
bar(1, nanmean(avgDist_prb_xbg(:,1) - avgDist_pre(:,1)), 'm');
bar(2, nanmean(avgDist_prb_xbg(:,2) - avgDist_pre(:,2)), 'r');
bar(3, nanmean(avgDist_prb_xbg(:,3) - avgDist_pre(:,3)), 'b');
bar(4, nanmean(avgDist_prb_xbg(:,4) - avgDist_pre(:,4)), 'k');
errorbar(1:4, nanmean(avgDist_prb_xbg - avgDist_pre), ...
    sqrt(nanvar(avgDist_prb_xbg - avgDist_pre, 0, 1)./sum(~isnan(avgDist_prb_xbg - avgDist_pre))), 'k.')
title('Probe Xbegin - Pre: diff');
axis([0 5 -.25 .05])

[~, p_prb_xbg, ~, t_prb_xbg] = ttest(avgDist_prb_xbg - avgDist_pre);
[~, p_prb_bg, ~, t_prb_bg] = ttest(avgDist_prb_bg - avgDist_pre);

[xan_1, xan_2, xan_3] = anova1(avgDist_prb_xbg - avgDist_pre);
[xc,xm,xh,xmns] = multcompare(xan_3);

[ban_1, ban_2, ban_3] = anova1(avgDist_prb_bg - avgDist_pre);
[bc,bm,bh,bmns] = multcompare(ban_3);
%%

figure;
plot([0 1], [0 1], '-', 'Color', [.5 .5 .5]);
hold on
plot(avgDist_pre(:, 1), avgDist_prb(:, 1), 'm.');
plot(avgDist_pre(:, 2), avgDist_prb(:, 2), 'r.');
plot(avgDist_pre(:, 3), avgDist_prb(:, 3), 'b.');
plot(avgDist_pre(:, 4), avgDist_prb(:, 4), 'k.');
lm1 = fitlm(avgDist_pre(:, 1), avgDist_prb(:, 1));
lm2 = fitlm(avgDist_pre(:, 2), avgDist_prb(:, 2));
lm3 = fitlm(avgDist_pre(:, 3), avgDist_prb(:, 3));
lm4 = fitlm(avgDist_pre(:, 4), avgDist_prb(:, 4));

plot([min(avgDist_pre(:,1)), max(avgDist_pre(:,1))], ...
    lm1.Coefficients.Estimate(1) + lm1.Coefficients.Estimate(2)*[min(avgDist_pre(:,1)),...
    max(avgDist_pre(:,1))], 'm-', 'LineWidth', 2.5);
plot([min(avgDist_pre(:,2)), max(avgDist_pre(:,2))], ...
    lm2.Coefficients.Estimate(1) + lm2.Coefficients.Estimate(2)*[min(avgDist_pre(:,2)),...
    max(avgDist_pre(:,2))], 'r-', 'LineWidth', 2.5);
plot([min(avgDist_pre(:,3)), max(avgDist_pre(:,3))], ...
    lm3.Coefficients.Estimate(1) + lm3.Coefficients.Estimate(2)*[min(avgDist_pre(:,3)),...
    max(avgDist_pre(:,3))], 'b-', 'LineWidth', 2.5);
plot([min(avgDist_pre(:,4)), max(avgDist_pre(:,4))], ...
    lm4.Coefficients.Estimate(1) + lm4.Coefficients.Estimate(2)*[min(avgDist_pre(:,4)),...
    max(avgDist_pre(:,4))], 'k-', 'LineWidth', 2.5);
axis([.3 1 .3 1])
% figure;
% hold on;
% plot([0 1], [0 0], 'k')
% plot(avgDist_pre(:, 1), lm1.Residuals.Raw, 'm.');
% plot(avgDist_pre(:, 2), lm2.Residuals.Raw, 'r.');
% plot(avgDist_pre(:, 3), lm3.Residuals.Raw, 'b.');
% plot(avgDist_pre(:, 4), lm4.Residuals.Raw, 'k.');

% %%
% 
% probes_binned{2} = mafilt(exact_track_dist_fall{2}(probe_trials{2}, :), 5);
% probes_binned{3} = mafilt(exact_track_dist_fall{3}(probe_trials{3}, :), 5);
% probes_binned{4} = mafilt(exact_track_dist_fall{4}(probe_trials{4}, :), 5);
% 
% probes_binned_all = [probes_binned{2}, probes_binned{3}, probes_binned{4}];
% probes_binned_diff = diff(probes_binned_all);
% [h_junk,p_bin_diff,c_bin_diff,stat_bin_diff] = ttest(probes_binned_diff');

%% show late and early performance across days:

% pre_inds_3 = {151:200, 351:400, [], 751:800};
pre_inds_3 = {151:200, 351:400, 551:600, 751:800};
post_inds_3 = {201:400, [401:475, 526:600], 601:800, 801:1000};
pre_inds_5 = {151:200, 351:400, 551:600, 751:800};
post_inds_5 = {201:400, 401:600, 601:800, [801:875, 926:1000]};
pre_inds_10 = {151:200, 351:400, 551:600, 751:800, ...
    951:1000, 1151:1200, 1351:1400, 1551:1600, 1751:1800};
post_inds_10 = {201:400, 401:600, 601:800, 801:1000, ...
    1001:1200, 1201:1400, 1401:1600, 1601:1800, [1801:1875, 1926:2000]};

pre_raw = nan(50, 20*3*2+20*5);
post_raw = nan(200, 20*3*2+20*5);

pre_buff = nan(50, 20*3*2+20*5);
post_buff = nan(200, 20*3*2+20*5);

figure;
for i_day = 1:2
    subplot(2,5,i_day); hold on;
    pre_ = [exact_track_dist_fall{2}(pre_inds_3{i_day}, :), ...
        exact_track_dist_fall{3}(pre_inds_5{i_day}, :), ...
        exact_track_dist_fall{4}(pre_inds_10{i_day}, :)];
    
    post_ = [[exact_track_dist_fall{2}(post_inds_3{i_day}, :); nan(200-length(post_inds_3{i_day}), size(exact_track_dist_fall{2},2))], ...
        [exact_track_dist_fall{3}(post_inds_5{i_day}, :); nan(200-length(post_inds_5{i_day}), size(exact_track_dist_fall{3},2))], ...
        [exact_track_dist_fall{4}(post_inds_10{i_day}, :); nan(200-length(post_inds_10{i_day}), size(exact_track_dist_fall{4},2))]];
    
    sub_norm_pre = nanmean(pre_./repmat(nanmean(pre_,1), size(pre_,1), 1), 2);
    plot((1:50)-51, sub_norm_pre);
    sub_norm_post = nanmean(post_./repmat(nanmean(pre_,1), size(post_,1),1),2);
    plot((51:250)-51, sub_norm_post, 'r');
    plot(([50 51])-51, [sub_norm_pre(end), sub_norm_post(1)], '--', 'Color', [.5 .5 .5]);
    axis([-50 200, .58 1.3])
end
i_pre = 1; i_post = 1;
for i_day = 3:4
    subplot(2,5,i_day); hold on;
    pre_ = [exact_track_dist_fall{2}(pre_inds_3{i_day}, :), ...
        exact_track_dist_fall{3}(pre_inds_5{i_day}, :), ...
        exact_track_dist_fall{4}(pre_inds_10{i_day}, :)];
    
    temp_3 = exact_track_dist_fall{2}(post_inds_3{i_day}, :);
    temp_5 = exact_track_dist_fall{3}(post_inds_5{i_day}, :);
    
    post_ = [[temp_3(1:(length(temp_3)/2),:); nan(200-length(post_inds_3{i_day}), size(exact_track_dist_fall{2},2)); temp_3((length(temp_3)/2+1):end, :)], ...
        [temp_5(1:(length(temp_5)/2),:); nan(200-length(post_inds_5{i_day}), size(exact_track_dist_fall{3},2)); temp_5((length(temp_5)/2+1):end, :)], ...
        [exact_track_dist_fall{4}(post_inds_10{i_day}, :); nan(200-length(post_inds_10{i_day}), size(exact_track_dist_fall{4},2))]];
    
    pre_raw(:, i_pre:(i_pre + size(pre_,2) - 1)) = pre_;
    post_raw(:, i_post:(i_post + size(post_,2) - 1)) = post_;
    
    pre_buff(:, i_pre:(i_pre + size(pre_,2) - 1)) = pre_./repmat(nanmean(pre_,1), size(pre_,1), 1);
    post_buff(:, i_post:(i_post + size(post_,2) - 1)) = post_./repmat(nanmean(pre_,1), size(post_,1), 1);
    
%     plot((1:50)-51, nanmean(pre_,2)./(nanmean(pre_,1)));
%     plot((51:250)-51, nanmean(post_,2)./(nanmean(pre_,1)), 'r');
%     plot(([50 51])-51, [nanmean(pre_(50, :),2), nanmean(post_(1, :),2)]./(nanmean(pre_,1)), '--', 'Color', [.5 .5 .5]);
    sub_norm_pre = nanmean(pre_./repmat(nanmean(pre_,1), size(pre_,1), 1), 2);
    plot((1:50)-51, sub_norm_pre);
    sub_norm_post = nanmean(post_./repmat(nanmean(pre_,1), size(post_,1),1),2);
    plot((51:250)-51, sub_norm_post, 'r');
    plot(([50 51])-51, [sub_norm_pre(end), sub_norm_post(1)], '--', 'Color', [.5 .5 .5]);
    axis([-50 200, .58 1.3])
    
    i_pre = i_pre + size(pre_,2);
    i_post = i_post + size(post_,2);
end
for i_day = 5:9
    subplot(2,5,i_day); hold on;
    pre_ = exact_track_dist_fall{4}(pre_inds_10{i_day}, :);
    
    temp_10 = exact_track_dist_fall{4}(post_inds_10{i_day}, :);
    
    post_ = [temp_10(1:(length(temp_10)/2),:); nan(200-length(post_inds_10{i_day}), size(exact_track_dist_fall{4},2)); temp_10((length(temp_10)/2+1):end, :)];
    
    pre_raw(1:size(pre_,1), i_pre:(i_pre + size(pre_,2) - 1)) = pre_;
    post_raw(1:size(post_,1), i_post:(i_post + size(post_,2) - 1)) = post_;
    
    pre_buff(1:size(pre_,1), i_pre:(i_pre + size(pre_,2) - 1)) = pre_./repmat(nanmean(pre_,1), size(pre_,1), 1);
    post_buff(1:size(post_,1), i_post:(i_post + size(post_,2) - 1)) = post_./repmat(nanmean(pre_,1), size(post_,1), 1);
    
%     plot((1:50)-51, nanmean(pre_,2)/nanmean(nanmean(pre_,1)));
%     plot((51:(50+size(post_,1)))-51, nanmean(post_,2)./(nanmean(pre_,1)), 'r');
%     plot(([50 51])-51, [nanmean(pre_(50, :),2), nanmean(post_(1, :),2)]./(nanmean(pre_,1)), '--', 'Color', [.5 .5 .5]);
    sub_norm_pre = nanmean(pre_./repmat(nanmean(pre_,1), size(pre_,1), 1), 2);
    plot((1:50)-51, sub_norm_pre);
    sub_norm_post = nanmean(post_./repmat(nanmean(pre_,1), size(post_,1),1),2);
    plot((51:(50+size(post_,1)))-51, sub_norm_post, 'r');
    plot(([50 51])-51, [sub_norm_pre(end), sub_norm_post(1)], '--', 'Color', [.5 .5 .5]);
    axis([-50 200, .58 1.3])
    
    i_pre = i_pre + size(pre_,2);
    i_post = i_post + size(post_,2);
end

%%
regress_window = 1:10;
early_window = 11:60;
late_window = 150:200;

figure; hold on;
plot((1:50)-51, nanmean(pre_buff,2));
plot((51:250)-51, nanmean(post_buff,2));
plot(([50 51])-51, [nanmean(pre_buff(end,:),2), nanmean(post_buff(1,:),2)], '--', 'Color', [.5 .5 .5]);
plot([-50 200], [1 1], '--b')
plot([150 200], repmat(nanmean(nanmean(post_buff(150:200, :),1)), 1,2), '--r')
plot([16, 65], repmat(nanmean(nanmean(post_buff(6:55, :),1),2), 1,2), '--m')

[h_junk, t_regress, c_regress, d_regress] = ttest(nanmean(post_buff(regress_window, :),1) - nanmean(pre_buff, 1));
[h_junk, t_overnight, c_over, d_over] = ttest(nanmean(post_buff(early_window, :),1) - nanmean(pre_buff, 1));
[h_junk, t_within, c_within, d_within] = ttest(nanmean(post_buff(late_window, :),1) - nanmean(post_buff(early_window, :),1));

% z_early = nan(150, size(post_raw,2));
% for i_sub = 1:size(post_raw,2)
%     for i_tr = 1:150
%         z_temp = zscore([post_raw(151:200, i_sub); post_raw(i_tr, i_sub)]);
%          z_early(i_tr, i_sub) = z_temp(end);
%     end
% end
% figure; hold on;
% plot(nanmean(z_early, 2));

%% determine largest learning rate cutoff
sum_bin_size = 4;
binned_raw = mafilt(post_raw, 5);
probes_diff = nan(size(binned_raw,1) - sum_bin_size, size(binned_raw,2));
for i_bin = 1:size(probes_diff,1)
    probes_diff(i_bin, :) = binned_raw(i_bin,:) - nanmean(binned_raw(i_bin+(1:4), :), 1);
end
% probes_binned_diff = diff(binned_raw);
[h_junk,p_bin_diff,c_bin_diff,stat_bin_diff2] = ttest(probes_diff');
p_new = pval_adjust(p_bin_diff, 'holm');

n_bins_different = 0; i_bin = 1;
sig_diff = p_new(i_bin) < 0.05;
while sig_diff
    sig_diff = p_new(i_bin) < 0.05;
    n_bins_different = n_bins_different + 1;
    i_bin = 1 + i_bin;
end
n_bins_different = n_bins_different - 1;

%%
rel_probe_signals = {nan(50, 21), nan(50, 20), nan(50, 20), nan(50, 20)};
asym_basis_signals = {nan(50, 21), nan(50, 20), nan(50, 20), nan(50, 20)};
probe_day_temp = [1, 3, 5, 10];
for i_group = 1:4
    %for i_sub = 1:size(exact_track_dist_fall{i_group},2)
        temp_asym = nanmean(exact_track_dist_fall{i_group}(preProbe_trials{i_group}, :), 1);
        rel_probe_signals{i_group} = exact_track_dist_fall{i_group}(probe_trials{i_group}, :)./repmat(temp_asym, length(probe_trials{i_group}), 1); 
        asym_basis_signals{i_group} = exact_track_dist_fall{i_group}(preProbe_trials{i_group}, :);
end

plt_clr = {'m', 'r', 'b', 'k'};
h1 = figure; hold on;
plot((1:5:50)-51, mafilt(nanmean([pre_buff, ...
    [asym_basis_signals{1}, asym_basis_signals{2}, asym_basis_signals{3}, asym_basis_signals{4}]./...
    repmat(nanmean([asym_basis_signals{1}, asym_basis_signals{2}, asym_basis_signals{3}, asym_basis_signals{4}],1), 50, 1)...
    ],2), 5));
plot((51:5:250)-51, mafilt(nanmean(post_buff,2),5));
for i_group = 1:4
    plot((51:5:100)-51, mafilt(nanmean(rel_probe_signals{i_group},2),5), plt_clr{i_group});
end
plot([-50 200], [1 1], '--b')
set(h1, 'Position', [440 378 675, 380])


%% make within-day and across-day learning curve demonstrations:
% early_temp = 1:50;
% late_temp = 151:200;
% % init_inds = {1:10, 200+(1:10), 400+(1:10), 600+(1:10), 800+(1:10), 1000+(1:10), 1200+(1:10), 1400+(1:10), 1600+(1:10), 1800+(1:10)};
% % early_inds = {11:60, 200+(11:60), 400+(11:60), 600+(11:60), 800+(11:60), 1000+(11:60), 1200+(11:60), 1400+(11:60), 1600+(11:60), 1800+(11:60)};
% early_inds = {early_temp, 200+(early_temp), 400+(early_temp), 600+(early_temp), 800+(early_temp), 1000+(early_temp), 1200+(early_temp), 1400+(early_temp), 1600+(early_temp), 1800+(early_temp)};
% late_inds = {late_temp, 200+(late_temp), 400+(late_temp), 600+(late_temp), 800+(late_temp), 1000+(late_temp), 1200+(late_temp), 1400+(late_temp), 1600+(late_temp), 1800+(late_temp)};
% ^- for window size of 50.

early_temp = 1:25;
late_temp = 176:200;
% init_inds = {1:10, 200+(1:10), 400+(1:10), 600+(1:10), 800+(1:10), 1000+(1:10), 1200+(1:10), 1400+(1:10), 1600+(1:10), 1800+(1:10)};
% early_inds = {11:60, 200+(11:60), 400+(11:60), 600+(11:60), 800+(11:60), 1000+(11:60), 1200+(11:60), 1400+(11:60), 1600+(11:60), 1800+(11:60)};
early_inds = {early_temp, 200+(early_temp), 400+(early_temp), 600+(early_temp), 800+(early_temp), 1000+(early_temp), 1200+(early_temp), 1400+(early_temp), 1600+(early_temp), 1800+(early_temp)};
late_inds = {late_temp, 200+(late_temp), 400+(late_temp), 600+(late_temp), 800+(late_temp), 1000+(late_temp), 1200+(late_temp), 1400+(late_temp), 1600+(late_temp), 1800+(late_temp)};
% ^- for window size of 25.

% init_day_buff = nan(10, 81, 10);
% late_day_buff = nan(50, 81, 10);
% early_day_buff = nan(50, 81, 10);
% ^- for window size of 50.
late_day_buff = nan(25, 81, 10);
early_day_buff = nan(25, 81, 10);
% ^- for window size of 25.

for i_day = 1:10
    ind_min = min(early_inds{i_day});
    k_sub = 1;
%     temp_day_init = nan(10, 81);
%     temp_day_early = nan(50, 81);
%     temp_day_late = nan(50, 81);
%       ^- for window size of 50.
    temp_day_early = nan(25, 81);
    temp_day_late = nan(25, 81);
    %   ^- for window size of 25.
    for i_group = 1:length(exact_track_dist_fall)
        if size(exact_track_dist_fall{i_group},1) > ind_min
%             temp_day_init(:, k_sub:(k_sub - 1 + size(exact_track_dist_fall{i_group},2))) = exact_track_dist_fall{i_group}(init_inds{i_day}, :);
            temp_day_early(:, k_sub:(k_sub - 1 + size(exact_track_dist_fall{i_group},2))) = exact_track_dist_fall{i_group}(early_inds{i_day}, :);
            temp_day_late(:, k_sub:(k_sub - 1 + size(exact_track_dist_fall{i_group},2))) = exact_track_dist_fall{i_group}(late_inds{i_day}, :);
            k_sub = k_sub + size(exact_track_dist_fall{i_group},2);
        end
    end
%     init_day_buff(:, 1:size(temp_day_init,2), i_day) = temp_day_init;
    early_day_buff(:, 1:size(temp_day_early,2), i_day) = temp_day_early;
    late_day_buff(:, 1:size(temp_day_late,2), i_day) = temp_day_late;
end

%% Test for change across days and within days:
% init_buff = reshape(nanmean(init_day_buff, 1), 81, 10);
early_buff = reshape(nanmean(early_day_buff, 1), 81, 10);
late_buff = reshape(nanmean(late_day_buff, 1), 81, 10);

% init_buff_3_10 = init_buff(22:end, 3:10);
early_buff_3_10 = early_buff(22:end, 3:10);
late_buff_3_10 = late_buff(22:end, 3:10);

diff_wDay = late_buff_3_10 - early_buff_3_10; %within day difference
diff_aDay = early_buff_3_10(:, 2:end) - late_buff_3_10(:, 1:(end-1)); %btwn day diff
% diff_decay = init_buff_3_10(:, 2:end) - late_buff_3_10(:, 1:(end-1));

diff_wDay_list = reshape(diff_wDay, 60*8, 1);
diff_aDay_list = reshape(diff_aDay, 60*7, 1);
% diff_dDay_list = reshape(diff_decay, 60*7, 1);

aDay_dat = repmat([(1:20)'; (21:40)'; (41:60)'], 7, 1);
wDay_dat = repmat([(1:20)'; (21:40)'; (41:60)'], 8, 1);
% dDay_dat = repmat([(1:20)'; (21:40)'; (41:60)'], 7, 1);

csvwrite('withinDay_design', wDay_dat);
csvwrite('withinDay_data', diff_wDay_list);
csvwrite('acrossDay_design', aDay_dat);
csvwrite('acrossDay_data', diff_aDay_list);
% csvwrite('decayDay_design', dDay_dat);
% csvwrite('decayDay_data', diff_dDay_list);

%%
figure;
% plot(0.05:9.05, reshape(nanmean(nanmean(init_day_buff,1),2), 10,1), 'r.')
hold on
plot(.3:9.3, reshape(nanmean(nanmean(early_day_buff,1),2), 10,1), 'k.')
plot((1:10), reshape(nanmean(nanmean(late_day_buff,1),2), 10,1), 'k.')
% plot([0.05:9.05; .3:9.3], ...
%     [reshape(nanmean(nanmean(init_day_buff,1),2), 1,10); reshape(nanmean(nanmean(early_day_buff,1),2), 1,10)],...
%     'k-');
% plot([0.3:9.3; 1:10], ...
%     [reshape(nanmean(nanmean(early_day_buff,1),2), 1,10); reshape(nanmean(nanmean(late_day_buff,1),2), 1,10)],...
%     'k-');
% plot([1:9; 1.05:9.05], ...
%     [reshape(nanmean(nanmean(late_day_buff(:,:,1:(end-1)),1),2), 1,9); ...
%     reshape(nanmean(nanmean(init_day_buff(:,:,2:end),1),2), 1,9)],...
%     'k:');
plot([0.3:9.3; 1:10], ...
    [reshape(nanmean(nanmean(early_day_buff,1),2), 1,10); reshape(nanmean(nanmean(late_day_buff,1),2), 1,10)],...
    'k-');
plot([1:9; 1.3:9.3], ...
    [reshape(nanmean(nanmean(late_day_buff(:,:,1:(end-1)),1),2), 1,9); ...
    reshape(nanmean(nanmean(early_day_buff(:,:,2:end),1),2), 1,9)],...
    'k:');

% init_day_mean = nan(81, 10);
early_day_mean = nan(81, 10);
late_day_mean = nan(81, 10);
for i_day = 1:size(late_day_buff,3)
%     temp_init_dat = init_day_buff(:,:,i_day);
    temp_early_dat = early_day_buff(:,:,i_day);
    temp_late_dat = late_day_buff(:,:,i_day);
    
%     errorbar(i_day -1 + 0.05, nanmean(nanmean(temp_init_dat,1),2), sqrt(nanvar(nanmean(temp_init_dat,1),0,2)/sum(~isnan(temp_init_dat(1,:)),2)), 'k.');
    errorbar(i_day -1 + .3, nanmean(nanmean(temp_early_dat,1),2), sqrt(nanvar(nanmean(temp_early_dat,1),0,2)/sum(~isnan(temp_early_dat(1,:)),2)), 'k.');
    errorbar(i_day, nanmean(nanmean(temp_late_dat,1),2), sqrt(nanvar(nanmean(temp_late_dat,1),0,2)/sum(~isnan(temp_late_dat(1,:)),2)), 'k.');
    
%     init_day_mean(:, i_day) = nanmean(temp_init_dat, 1)';
    early_day_mean(:, i_day) = nanmean(temp_early_dat, 1)';
    late_day_mean(:, i_day) = nanmean(temp_late_dat, 1)';
end

% assess within day learning:
pooled_early_late_dist = [early_day_mean(:); late_day_mean(:)];
pooled_early_late_category = [repmat((1:81)', 20, 1), [zeros(81*10, 1); ones(81*10, 1)],...
    repmat([ones(81, 1); 2*ones(81, 1); 3*ones(81, 1); 4*ones(81, 1); 5*ones(81, 1);...
    6*ones(81, 1); 7*ones(81, 1); 8*ones(81, 1); 9*ones(81, 1); 10*ones(81, 1)], 2, 1)];

csvwrite('pooled_early_late_dist', pooled_early_late_dist);
csvwrite('pooled_early_late_category', pooled_early_late_category);

T = table;
T.subject = pooled_early_late_category(:,1);
T.window = pooled_early_late_category(:,2);
T.day = pooled_early_late_category(:,3);
T.response = pooled_early_late_dist;
lme_win = fitlme(T, 'response ~ window + day + day:window + (1|subject)');
lme_win__x = fitlme(T, 'response ~ window + day + (1|subject)');
lme_win__w = fitlme(T, 'response ~ day + (1|subject)');
lme_win__d = fitlme(T, 'response ~ window + (1|subject)');

[av_win_x, sim_inf_x] = compare(lme_win, lme_win__x, 'nsim', 200);
[av_win_w, sim_inf_w] = compare(lme_win, lme_win__w, 'nsim', 200);
[av_win_d, sim_inf_d] = compare(lme_win, lme_win__d, 'nsim', 200);

anova_win = anova(lme_win);

% assess overnight learning:
pooled_early_late_dist = [reshape(late_day_mean(:, 1:(end-1)), 81*9, 1); reshape(early_day_mean(:, 2:end), 81*9, 1)];
pooled_early_late_category = [repmat((1:81)', 18, 1), [zeros(81*9, 1); ones(81*9, 1)],...
    repmat([ones(81, 1); 2*ones(81, 1); 3*ones(81, 1); 4*ones(81, 1); 5*ones(81, 1);...
    6*ones(81, 1); 7*ones(81, 1); 8*ones(81, 1); 9*ones(81, 1)], 2, 1)];

csvwrite('pooled_overnight_dist', pooled_early_late_dist);
csvwrite('pooled_overnight_category', pooled_early_late_category);

% assess decay across days:
% pooled_init_late_dist = [reshape(late_day_mean(:, 1:(end-1)), 81*9, 1); reshape(init_day_mean(:, 2:end), 81*9, 1)];
% pooled_init_late_category = [repmat((1:81)', 18, 1), [zeros(81*9, 1); ones(81*9, 1)],...
%     repmat([ones(81, 1); 2*ones(81, 1); 3*ones(81, 1); 4*ones(81, 1); 5*ones(81, 1);...
%     6*ones(81, 1); 7*ones(81, 1); 8*ones(81, 1); 9*ones(81, 1)], 2, 1)];

% csvwrite('pooled_overnight_dist', pooled_init_late_dist);
% csvwrite('pooled_overnight_category', pooled_init_late_category);


%%
within_day_diff = late_day_buff - early_day_buff;
btwn_day_diff = diff(late_day_buff, 1, 3);

figure; hold on
plot(1:10, reshape(nanmean(nanmean(within_day_diff,1),2),10,1), '.-');
plot(1.5:10, reshape(nanmean(nanmean(btwn_day_diff,1),2),9,1), '.-');
    