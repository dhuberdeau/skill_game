load('traj_pcData.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat')
load('exact_track_dist_full_v3.mat') %distance travelled measure
% load('score_slopes_data')
load('score_by_trial_grp')
%% collect trajectory data, distance travelled data, and policy deviation data 
% from each day's "window A" (pre- on probe days), and "window B" (probe-
% on probe days). Probe window data is in it's own structure.. "*_probe_*",
% while training data is in "*_winds_*"

pre_wind = [26:75; 226:275; 426:475; 626:675; 826:875; 1026:1075; 1226:1275; 1426:1475; 1626:1675; 1826:1875]';
probe_wind = [76:125; 276:325; 476:525; 676:725; 876:925; 1076:1125; 1276:1325; 1476:1525; 1676:1725; 1876:1925]';

subs_grp = {1:21, 22:41, 43:62, 64:83};
prb_wind_grp = [2, 6, 10, 20];
% prb_wind_grp = [1, 5, 9, 19];
n_wind_grp = [2, 10, 10, 20];
n_days_grp = [1, 5, 5, 10];
wind_len = 50;
n_sub = 20;

n_over_samp = [0, 0, 0, 0];
n_over_dt = [0, 0, 0, 0];
n_over_dev = [0, 0, 0, 0];

samp_winds_grp = cell(1,4);
samp_probe_grp = cell(1,4);
dev_winds_grp = cell(1,4);
dev_probe_grp = cell(1,4);
dt_winds_grp = cell(1,4);
dt_probe_grp = cell(1,4);
for i_grp = 2:4
    % reshape pcData:
    samp_winds = nan(wind_len, n_wind_grp(i_grp) - 1, 20); %trials X windows X subjects
    samp_probe = nan(wind_len, n_sub); %trials X subjects
    for i_sub = 1:n_sub
        temp_pre = reshape(traj_pcData_pre{1}(subs_grp{i_grp}(i_sub), 1:(n_wind_grp(i_grp)/2), :), n_wind_grp(i_grp)/2, wind_len)';
        temp_prb = reshape(traj_pcData_prb{1}(subs_grp{i_grp}(i_sub), 1:(n_wind_grp(i_grp)/2), :), n_wind_grp(i_grp)/2, wind_len)';
        
        temp_all = reshape([temp_pre; temp_prb], wind_len, n_wind_grp(i_grp));
        samp_winds(:, 1:(n_wind_grp(i_grp) - 1), i_sub) = temp_all(:, setdiff(1:n_wind_grp(i_grp), prb_wind_grp(i_grp)));
        samp_probe(:, i_sub) = temp_all(:, prb_wind_grp(i_grp));
    end
    
    % reshape dist travel data:
    dt_winds = nan(wind_len, n_wind_grp(i_grp) - 1, 20);
    dt_probe = nan(wind_len, n_sub);
    for i_sub = 1:n_sub
        temp_k = pre_wind(:, 1:n_days_grp(i_grp));
        temp_pre = reshape(exact_track_dist_fall{i_grp}(temp_k(:), i_sub), wind_len, n_wind_grp(i_grp)/2);
        temp_k = probe_wind(:, 1:n_days_grp(i_grp));
        temp_prb = reshape(exact_track_dist_fall{i_grp}(temp_k(:), i_sub), wind_len, n_wind_grp(i_grp)/2);
        
        temp_all = reshape([temp_pre; temp_prb], wind_len, n_wind_grp(i_grp));
        dt_winds(:, 1:(n_wind_grp(i_grp) - 1), i_sub) = temp_all(:, setdiff(1:n_wind_grp(i_grp), prb_wind_grp(i_grp)));
        dt_probe(:, i_sub) = temp_all(:, prb_wind_grp(i_grp));
    end
    
    % reshape deviation data:
    dev_winds = nan(wind_len, n_wind_grp(i_grp) - 1, 20);
    dev_probe = nan(wind_len, n_sub);
    for i_sub = 1:n_sub
        temp_pre = reshape(score_succ_by_trial_grp{1, i_grp}(:, i_sub, :), wind_len, n_wind_grp(i_grp)/2);
        temp_prb = reshape(score_succ_by_trial_grp{2, i_grp}(:, i_sub, :), wind_len, n_wind_grp(i_grp)/2);
        
        temp_all = reshape([temp_pre; temp_prb], wind_len, n_wind_grp(i_grp));
        dev_winds(:, 1:(n_wind_grp(i_grp) - 1), i_sub) = temp_all(:, setdiff(1:n_wind_grp(i_grp), prb_wind_grp(i_grp)));
        dev_probe(:, i_sub) = temp_all(:, prb_wind_grp(i_grp));
    end
    
    
    samp_winds_grp{i_grp} = samp_winds;
    samp_probe_grp{i_grp} = samp_probe;
    dt_winds_grp{i_grp} = dt_winds;
    dt_probe_grp{i_grp} = dt_probe;
    dev_winds_grp{i_grp} = dev_winds;
    dev_probe_grp{i_grp} = dev_probe;
end
 
%% 
ks_samp_grp = cell(1,4);
ks_dt_grp = cell(1,4);
ks_dev_grp = cell(1,4);
match_samp_grp = cell(1,4);
match_dt_grp = cell(1,4);
match_dev_grp = cell(1,4);
for i_grp = 2:4

    if i_grp == 4
        h1_a = figure;
        h2_a = figure;
        h1_b = figure;
        h2_b = figure;
        h1_c = figure;
        h2_c = figure;
    end
    
    temp_samp = nan(n_sub,1);
    temp_dt = nan(n_sub,1);
    temp_dev = nan(n_sub,1);
    for i_sub = 1:n_sub
        x_temp = 1:(n_wind_grp(i_grp) - 1);
        
        try
        m_temp = nanmean(samp_winds_grp{i_grp}(:, :, i_sub), 1);
        s_temp = nanstd(samp_winds_grp{i_grp}(:, :, i_sub), 0, 1);
        lm_m = fitlm(x_temp, m_temp);
        a_m = lm_m.Coefficients.Estimate(1); b_m = lm_m.Coefficients.Estimate(2);
        lm_s = fitlm(x_temp, s_temp);
        a_s = lm_s.Coefficients.Estimate(1); b_s = lm_s.Coefficients.Estimate(2);
        if i_grp == 4
            figure(h1_a);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), m_temp, '.');
            plot(x_temp(1:19), a_m + b_m*x_temp(1:19), 'k-');
            title('m:samp')
            figure(h2_a);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), s_temp, '.')
            plot(x_temp(1:19), a_s + b_s*x_temp(1:19), 'k-');
            title('S:samp')
        end
        f = @(x_i,w) (1./sqrt(2*pi*(a_s + b_s*w).^2)).*exp(-(x_i - (a_m + b_m*w)).^2/(2*(a_s + b_s*w).^2));
        temp_samp(i_sub) = mle(samp_probe_grp{i_grp}(~isnan(samp_probe_grp{i_grp}(:, i_sub)), i_sub), 'pdf', f, 'start', 20);
        catch end
        
        try
        m_temp = nanmean(dt_winds_grp{i_grp}(:, :, i_sub), 1);
        s_temp = nanstd(dt_winds_grp{i_grp}(:, :, i_sub), 0, 1);
        lm_m = fitlm(x_temp, m_temp);
        a_m = lm_m.Coefficients.Estimate(1); b_m = lm_m.Coefficients.Estimate(2);
        lm_s = fitlm(x_temp, s_temp);
        a_s = lm_s.Coefficients.Estimate(1); b_s = lm_s.Coefficients.Estimate(2);
        if i_grp == 4
            figure(h1_b);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), m_temp, '.');
            plot(x_temp(1:19), a_m + b_m*x_temp(1:19), 'k-');
            title('m:dt')
            figure(h2_b);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), s_temp, '.')
            plot(x_temp(1:19), a_s + b_s*x_temp(1:19), 'k-');
            title('S:dt')
        end
        f = @(x_i,w) (1./sqrt(2*pi*(a_s + b_s*w).^2)).*exp(-(x_i - (a_m + b_m*w)).^2/(2*(a_s + b_s*w).^2));
        temp_dt(i_sub) = mle(dt_probe_grp{i_grp}(~isnan(dt_probe_grp{i_grp}(:, i_sub)), i_sub), 'pdf', f, 'start', 20);
        catch end
        
        try
        m_temp = nanmean(dev_winds_grp{i_grp}(:, :, i_sub), 1);
        s_temp = nanstd(dev_winds_grp{i_grp}(:, :, i_sub), 0, 1);
        lm_m = fitlm(x_temp, m_temp);
        a_m = lm_m.Coefficients.Estimate(1); b_m = lm_m.Coefficients.Estimate(2);
        lm_s = fitlm(x_temp, s_temp);
        a_s = lm_s.Coefficients.Estimate(1); b_s = lm_s.Coefficients.Estimate(2);
        if i_grp == 4
            figure(h1_c);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), m_temp, '.');
            plot(x_temp(1:19), a_m + b_m*x_temp(1:19), 'k-');
            title('m:dev')
            figure(h2_c);
            subplot(4,5,i_sub); hold on;
            plot(x_temp(1:19), s_temp, '.')
            plot(x_temp(1:19), a_s + b_s*x_temp(1:19), 'k-');
            title('S:dev')
        end
        f = @(x_i,w) (1./sqrt(2*pi*(a_s + b_s*w).^2)).*exp(-(x_i - (a_m + b_m*w)).^2/(2*(a_s + b_s*w).^2));
        temp_dev(i_sub) = mle(dev_probe_grp{i_grp}(~isnan(dev_probe_grp{i_grp}(:, i_sub)), i_sub), 'pdf', f, 'start', 20);
        catch end
    end
    
    temp_samp(temp_samp == n_wind_grp(i_grp)) = nan;
    temp_samp(temp_samp < 0) = 0;
    n_over_samp(i_grp) = n_over_samp(i_grp) + sum(temp_samp > n_wind_grp(i_grp));
    temp_samp(temp_samp > n_wind_grp(i_grp)) = n_wind_grp(i_grp);
    
    temp_dt(temp_dt == n_wind_grp(i_grp)) = nan;
    temp_dt(temp_dt < 0) = 0;
    n_over_dt(i_grp) = n_over_dt(i_grp) + sum(temp_dt > n_wind_grp(i_grp));
    temp_dt(temp_dt > n_wind_grp(i_grp)) = n_wind_grp(i_grp);
    
    temp_dev(temp_dev == n_wind_grp(i_grp)) = nan;
    temp_dev(temp_dev < 0) = 0;
    n_over_dev(i_grp) = n_over_dev(i_grp) + sum(temp_dev > n_wind_grp(i_grp));
    temp_dev(temp_dev > n_wind_grp(i_grp)) = n_wind_grp(i_grp);
    
    match_samp_grp{i_grp} = n_wind_grp(i_grp) - temp_samp;
    match_dt_grp{i_grp} = n_wind_grp(i_grp) - temp_dt;
    match_dev_grp{i_grp} = n_wind_grp(i_grp) - temp_dev;

end

%%
match_samp_all = [match_samp_grp{2}; match_samp_grp{3}; match_samp_grp{4}];
match_dt_all = [match_dt_grp{2}; match_dt_grp{3}; match_dt_grp{4}];
match_dev_all = [match_dev_grp{2}; match_dev_grp{3}; match_dev_grp{4}];

match_samp_temp = match_samp_all(match_samp_all <= 10);
match_dt_temp = match_dt_all(match_dt_all <= 10);
match_dev_temp = match_dev_all(match_dev_all <= 10);

l_samp = poissfit(match_samp_temp(~isnan(match_samp_temp)));
l_dt = poissfit(match_dt_temp(~isnan(match_dt_temp)));
l_dev = poissfit(match_dev_temp(~isnan(match_dev_temp)));
% l_samp = betafit(match_samp_temp(~isnan(match_samp_temp)));
% l_dt = betafit(match_dt_temp(~isnan(match_dt_temp)));
% l_dev = betafit(match_dev_temp(~isnan(match_dev_temp)));

x = 0:1:20;
c = [x, 21];

figure;
subplot(3,1,1); hold on;
[n] = histc(match_samp_all, c);
bar(-x+20, n(1:length(x))/sum(n(1:length(x))));
plot(-x+20, poisspdf(x, l_samp));
% plot(-x+20, betapdf(x, l_samp(1), l_samp(2)));
axis([0 21 0 .4])

subplot(3,1,2); hold on;
[n] = histc(match_dev_all, c);
bar(-x+20, n(1:length(x))/sum(n(1:length(x))));
plot(-x+20, poisspdf(x, l_dev));
% plot(-x+20, betapdf(x, l_dev(1), l_dev(2)));
axis([0 21 0 .4])

subplot(3,1,3); hold on;
[n] = histc(match_dt_all, c);
bar(-x+20, n(1:length(x))/sum(n(1:length(x))));
plot(-x+20, poisspdf(x, l_dt));
% plot(-x+20, betapdf(x, l_dt(1), l_dt(2)));
axis([0 21 0 .4])

figure; hold on;
plot(-x+20, poisspdf(x, l_samp));
plot(-x+20, poisspdf(x, l_dev));
plot(-x+20, poisspdf(x, l_dt));
% plot(-x+20, betapdf(x, l_samp(1), l_samp(2)));
% plot(-x+20, betapdf(x, l_dev(1), l_dev(2)));
% plot(-x+20, betapdf(x, l_dt(1), l_dt(2)));
legend('Kinematics', 'Deviation', 'Dist Trv')

[h_samp_dt, p_samp_dt] = kstest2(match_samp_all, match_dt_all);
[h_samp_dev, p_samp_dev] = kstest2(match_samp_all, match_dev_all);
[h_dev_dt, p_dev_dt] = kstest2(match_dev_all, match_dt_all);

[h_anova, p_anova, stat_anova] = anova1([match_samp_all, match_dev_all, match_dt_all, ]);

kruk_p = kruskalwallis([match_samp_all, match_dev_all, match_dt_all, ])