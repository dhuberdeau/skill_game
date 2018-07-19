

load('tilt_dir.mat')
load('tilt_mag.mat')
load('traj_dir.mat')
load('traj_loc.mat')
load('exact_track_dist_full_v3.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat') % replaced with
% above, 9/27/17

fall_off_succ_TH = .8;
sample_size_th = 10;
z_scr_TH = 50;

if ~exist('i_grp', 'var')
    i_grp = 4;
end
grp_remove_sub = {[], [6], [2,3], [11 12 19]};
% REMOVE_SUB = [2 4 11 19 20];
REMOVE_SUB = grp_remove_sub{i_grp};
% REMOVE_SUB = [];

% loop_discard_TH = 200;
% temp_base = rad2deg(unwrap(deg2rad(temp_base)));
% traj_basis{i_day} = temp_base(:, sum(temp_base > loop_discard_TH, 1) < 1);

%%
% label successful trajectories all through segment
traj_succ = cell(size(traj_loc));

for k_grp = 1:length(traj_succ)
    traj_succ{k_grp} = nan(size(traj_loc{k_grp}));
    for i_sub = 1:size(traj_succ{k_grp},3)
        traj_succ{k_grp}(:,:,i_sub) = repmat(exact_track_dist_fall{k_grp}(:,i_sub)' > fall_off_succ_TH, size(traj_loc{k_grp},1), 1);
    end
end

%% compute Pr(visit location | succ)
day_inds = {...
    1:200, ...
    1:1000, ...
    1:1000, ...
    1:2000};
% day_inds = {[51:150], [251:350], [451:550], [651:750], [851:950],...
%     1000+[51:150], 1000+[251:350], 1000+[451:550], 1000+[651:750], 1000+[851:950]};

grp_days = [1 5 5 10];
n_subs = 20;

p_loc_g_succ = cell(1,n_subs);
u_loc = cell(2,n_subs);

num_l_bins = 10;
L_BINS = linspace(0,1,num_l_bins+1);
num_d_bins = 20; %number of direction bins
D_BINS = linspace(-180, 180, num_d_bins + 1);
num_locs = 20;
% NOTE: TRACK LENGTHS USED ARE .25 TO .45

% for i_day = 1:grp_days(i_grp)
sub_inds = 1:size(traj_succ{i_grp},3);
% moved * from here

n = nan(num_l_bins,num_d_bins,num_locs);
u = nan(num_l_bins,num_d_bins,num_locs);
u_s = nan(num_l_bins,num_d_bins,num_locs);

for i_sub = 1:n_subs
    % * moved to here
    succ_day_ = permute(traj_succ{i_grp}(:,day_inds{i_grp},setdiff(sub_inds, i_sub)), [2, 3, 1]);
    succ_day = reshape(succ_day_, size(succ_day_, 1)*size(succ_day_,2), size(succ_day_,3));

    loc_day_ = permute(traj_loc{i_grp}(:,day_inds{i_grp},setdiff(sub_inds, i_sub)), [2, 3, 1]);
    loc_day = reshape(loc_day_, size(loc_day_, 1)*size(loc_day_,2), size(succ_day_,3));

    tilt_day_ = permute(tilt_dir{i_grp}(:,day_inds{i_grp},setdiff(sub_inds, i_sub)), [2 3 1]);
    tilt_day = reshape(tilt_day_, size(tilt_day_, 1)*size(tilt_day_,2), size(succ_day_,3));

    dir_day_ = permute(traj_dir{i_grp}(:,day_inds{i_grp},setdiff(sub_inds, i_sub)), [2 3 1]);
    dir_day = reshape(dir_day_, size(dir_day_, 1)*size(dir_day_,2), size(succ_day_,3));
    try
        temp_loc = loc_day(succ_day(:,1)>0,:);
        temp_tilt = tilt_day(succ_day(:,1)>0,:);
        temp_dir = dir_day(succ_day(:,1)>0,:);
        
        temp_n = nan(num_l_bins, num_d_bins, num_locs);
        for i_loc = 1:size(temp_loc,2)
            temp_n_ = hist3([temp_loc(:,i_loc), temp_dir(:, i_loc)], {L_BINS, D_BINS});
            temp_n(:,:,i_loc) = temp_n_(1:num_l_bins,1:num_d_bins);
        end
        n(:, :, :, i_sub) = temp_n;
        for i_loc = 1:num_locs
            for i_bin = 1:(num_l_bins-1)
                for i_dir = 1:(num_d_bins-1)
                    crit1 = temp_loc(:,i_loc) > L_BINS(i_bin) & temp_loc(:,i_loc) <= L_BINS(i_bin + 1);
                    crit2 = temp_dir(:,i_loc) > D_BINS(i_dir) & temp_dir(:,i_loc) <= D_BINS(i_dir + 1);
                    if sum(crit1&crit2) >= sample_size_th
                        u(i_bin, i_dir, i_loc) = mean(temp_tilt(crit1 & crit2, i_loc));
                        u_s(i_bin, i_dir, i_loc) = var(temp_tilt(crit1 & crit2, i_loc));
                    end
                end
            end
        end
    catch
        er = lasterror;
        warning('wtf')
    end
    u_loc{1, i_sub} = u;
    u_loc{2, i_sub} = u_s;
end

disp(['Policy map complete']);
% end

%%
u_temp = nan([size(u_loc{1,1}), n_subs]);
for i_sub = 1:n_subs
    u_temp(:,:,:, i_sub) = u_loc{1, i_sub};
end
figure; hold on;
for i_bin = 1:size(u_temp,1)
    for i_dir = 1:size(u_temp,2)
        for i_loc = 1:size(u_temp,3)
            plot3(i_loc + .5*[0, cosd(nanmean(u_temp(i_bin, i_dir, i_loc, :),4))],...
                i_bin/size(u_temp,1) + .5*[0, sind(nanmean(u_temp(i_bin, i_dir, i_loc, :),4))],...
                [i_dir/size(u_temp,2), i_dir/size(u_temp,2)], 'k.-')
        end
    end
end

%% compute Pr(visit location | succ) for pre- probe window
% pre_wind = {1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050, 1201:1250, 1401:1450, 1601:1650, 1801:1850};
pre_wind = {1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050, 1201:1250, 1401:1450, 1601:1650, 1801:1850};

zu_trial = nan(length(pre_wind{1}), num_locs, n_subs, 10); %10 days

for i_day = 1:grp_days(i_grp)
    succ_day = permute(traj_succ{i_grp}(:,pre_wind{i_day},:), [2, 1, 3]);
    loc_day = permute(traj_loc{i_grp}(:,pre_wind{i_day},:), [2, 1, 3]);
    tilt_day = permute(tilt_dir{i_grp}(:,pre_wind{i_day},:), [2 1 3]);
    dir_day = permute(traj_dir{i_grp}(:,pre_wind{i_day},:), [2 1 3]);
    
    zu_day = nan(size(tilt_day));
%     zus_day = nan(num_l_bins, num_d_bins, num_locs, n_subs);
    for i_sub = 1:size(succ_day,3)
        trial_inds = 1:size(succ_day,1);
        target_trials = trial_inds(succ_day(:,1,i_sub)<1);
        
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
%% compute Pr(visit location | succ) for probe window
probe_wind = {151:200, 351:400, 551:600, 751:800, 951:1000, 1151:1200, 1351:1400, 1551:1600, 1751:1800, 1951:2000};
% probe_wind = {76:125, 276:325, 476:525, 676:725, 876:925, 1076:1125, 1276:1325, 1476:1525, 1676:1725, 1876:1925};

zu_trial_probe = nan(length(pre_wind{1}), num_locs, n_subs, grp_days(i_grp)); %10 days

for i_day = 1:grp_days(i_grp)
    succ_day = permute(traj_succ{i_grp}(:,probe_wind{i_day},:), [2, 1, 3]);
    loc_day = permute(traj_loc{i_grp}(:,probe_wind{i_day},:), [2, 1, 3]);
    tilt_day = permute(tilt_dir{i_grp}(:,probe_wind{i_day},:), [2 1 3]);
    dir_day = permute(traj_dir{i_grp}(:,probe_wind{i_day},:), [2 1 3]);
    
    zu_day = nan(size(tilt_day));
%     zus_day = nan(num_l_bins, num_d_bins, num_locs, n_subs);
    for i_sub = 1:size(succ_day,3)
        trial_inds = 1:size(succ_day,1);
        target_trials = trial_inds(succ_day(:,1,i_sub)<1);
        
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
                                mu = u_loc{1,i_sub}(i_bin, i_dir, i_loc, :);
%                                 s_temp = reshape(u_conglom{2,i_sub}(i_bin, i_dir, i_loc, :), 10, 1);
                                sig = sqrt(u_loc{2,i_sub}(i_bin, i_dir, i_loc, :));
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

%% 
zu_rev = nan(size(zu_trial));
zu_rev_probe = nan(size(zu_trial_probe));
loc_inds = 1:num_locs;
for i_day = 1:grp_days(i_grp)
    for i_sub = setdiff(1:size(zu_trial,3), REMOVE_SUB)
        for i_tr = 1:size(zu_trial,1)
            temp = zu_trial(i_tr, :, i_sub, i_day);
            k_max = max(loc_inds(~isnan(temp)));
            temp2 = temp(1:k_max);
            zu_rev(i_tr, (num_locs - k_max + 1):num_locs, i_sub, i_day) = temp2;
            
            temp = zu_trial_probe(i_tr, :, i_sub, i_day);
            k_max = max(loc_inds(~isnan(temp)));
            temp2 = temp(1:k_max);
            zu_rev_probe(i_tr, (num_locs - k_max + 1):num_locs, i_sub, i_day) = temp2;
        end
    end
end

figure;
for i_day = 1:grp_days(i_grp)
    subplot(2,5,i_day); hold on;
    dev_pre = nan(n_subs, num_locs);
    dev_prb = nan(n_subs, num_locs);
    for i_sub = setdiff(1:size(zu_rev,3), REMOVE_SUB)
        temp = zu_rev(:,:,i_sub, i_day);
        temp_pre = temp(sum(isnan(temp),2) < size(temp,2), :);

        temp = zu_rev_probe(:,:,i_sub, i_day);
        temp_prb = temp(sum(isnan(temp),2) < size(temp,2), :);
        
        dev_pre(i_sub, :) = nanmean(temp_pre,1);
        dev_prb(i_sub, :) = nanmean(temp_prb,1);
    end
    errorfield(1:num_locs, nanmean(dev_pre,1), nanstd(dev_pre)./sqrt(sum(~isnan(dev_pre(:,end)))), 'k');
    errorfield(1:num_locs, nanmean(dev_prb,1), nanstd(dev_prb)./sqrt(sum(~isnan(dev_prb(:,end)))), 'r');
    axis([0 20 0 5])
end
%%
pre_score = nan(n_subs, grp_days(i_grp));
prb_score = nan(n_subs, grp_days(i_grp));
for i_day = 1:grp_days(i_grp)
    for i_sub = setdiff(1:20, REMOVE_SUB)
        temp = nanmean(zu_rev(:, :, i_sub, i_day), 2);
        pre_score(i_sub, i_day) = nanmean(temp,1);

        temp_p = nanmean(zu_rev_probe(:, :, i_sub, i_day), 2);
        prb_score(i_sub, i_day) = nanmean(temp_p,1);
    end
end

pre_wind_wind = [[.375 1.375, 2.375, 3.375 4.375], 5+[.375 1.375, 2.375, 3.375 4.375]];
prb_wind_wind = [[.625 1.625, 2.625, 3.625 4.625], 5+[.625 1.625, 2.625, 3.625 4.625]];
x = reshape([pre_wind_wind(1:grp_days(i_grp)); prb_wind_wind(1:grp_days(i_grp))], 2*grp_days(i_grp), 1);

score = reshape([pre_score; prb_score], n_subs, 2*grp_days(i_grp));
% score(score > 6) = nan;


score(REMOVE_SUB, :) = nan;

figure; hold on;
errorbar(x, nanmean(score,1), nanstd(score)./sqrt(sum(~isnan(score),1)), 'k.-');

wind_probe = [nan, 6, 10, 20];
day_probe = [nan, 3, 5, 10];

[h_dev,p_dev] = ttest(score(:, wind_probe(i_grp) - 1) - score(:, wind_probe(i_grp)));

temp_diff = diff(score,1,2);
score_diff = temp_diff(:,1:2:(2*grp_days(i_grp) - 1));
[a_anv, b_anv, c_anv] = anova1(score_diff(:, [day_probe(i_grp)-2, day_probe(i_grp)-1, day_probe(i_grp)]));
figure; hold on;
plot([1 2], [score(:, wind_probe(i_grp) - 1), score(:, wind_probe(i_grp))], 'Color', [.8 .8 .8]);
errorbar([1 2], nanmean([score(:, wind_probe(i_grp)-1), score(:, wind_probe(i_grp))],1), nanstd([score(:, wind_probe(i_grp)-1), score(:, wind_probe(i_grp))])./sqrt(sum(~isnan([score(:, wind_probe(i_grp)-1), score(:, wind_probe(i_grp))]),1)), 'k.-');

%% fit linear model and test for effect of probe

sub_mat = repmat((1:20)', 1, 2*grp_days(i_grp));
win_mat = repmat(reshape([pre_wind_wind(1:grp_days(i_grp)); prb_wind_wind(1:grp_days(i_grp))], 2*grp_days(i_grp), 1)', 20, 1);
if i_grp == 2
    prb_mat = logical([zeros(20, 5), ones(20,1), zeros(20, 4)]);
else
    prb_mat = logical([zeros(20, 2*grp_days(i_grp) - 1), ones(20,1)]);
end

lab_mat = [reshape(sub_mat, 20*(2*grp_days(i_grp)), 1), reshape(win_mat, 20*(2*grp_days(i_grp)), 1), reshape(prb_mat, 20*(2*grp_days(i_grp)), 1)];
res_mat = reshape(score, 20*(2*grp_days(i_grp)), 1);
    
T = table;
T.subject = categorical(lab_mat(:, 1));
T.window = double(lab_mat(:, 2));
T.probe = logical(lab_mat(:, 3));
T.response = double(res_mat);

lm = fitlme(T, 'response ~ window + probe + (1|subject)');

%%
dt_pre = nan(20, grp_days(i_grp));
dt_prb = nan(20, grp_days(i_grp));

for i_day = 1:grp_days(i_grp)
    dt_pre(:, i_day) = nanmean(exact_track_dist_fall{i_grp}(pre_wind{i_day}, :), 1)';
    dt_prb(:, i_day) = nanmean(exact_track_dist_fall{i_grp}(probe_wind{i_day}, :), 1)';
end

dist_mat = reshape([dt_pre; dt_prb], 20, 2*grp_days(i_grp));
T.dist = double(reshape(dist_mat, 20*(2*grp_days(i_grp)), 1));

x_fit = log(abs(T.response(T.window < (grp_days(i_grp)-1))));
x_fit(x_fit < -3.5 | x_fit > 3.5) = nan;
y_fit = T.dist(T.window < (grp_days(i_grp)-1));
y_fit(isnan(x_fit)) = nan;

pre_wind_grp = [nan, 2.3750, 4.3750, 9.3750];
probe_wind_grp = [nan, 2.6250, 4.6250, 9.6250];

figure; hold on;
plot(x_fit, y_fit, '.', 'Color', [.7 .7 .7]);
plot(log(abs(T.response(T.window == pre_wind_grp(i_grp) & ~T.probe))), T.dist(T.window == pre_wind_grp(i_grp) & ~T.probe), 'k.', 'MarkerSize', 18);
plot(log(abs(T.response(T.window == probe_wind_grp(i_grp) & T.probe))), T.dist(T.window == probe_wind_grp(i_grp) & T.probe), 'r.', 'MarkerSize', 18);

% fitmod = fitlm(log(abs(T.response(T.window < 9))), T.dist(T.window < 9), 'RobustOpts', 'on');
fitmod = fitlm(x_fit, y_fit, 'RobustOpts', 'on');
fitfunc = @(t1) fitmod.Coefficients.Estimate(1) + fitmod.Coefficients.Estimate(2)*t1;

t = [min(x_fit), max(x_fit)];
l = fitfunc(t);
plot(t, l)

resid_pre = (fitfunc(log(abs(T.response(T.window == pre_wind_grp(i_grp) & ~T.probe)))) - T.dist(T.window == pre_wind_grp(i_grp) & ~T.probe));
resid_prb = (fitfunc(log(abs(T.response(T.window == probe_wind_grp(i_grp) & T.probe)))) - T.dist(T.window == probe_wind_grp(i_grp) & T.probe));

figure; hold on;
plot([1 2], [resid_pre, resid_prb], 'Color', [.8 .8 .8]);
errorbar([1 2], nanmean([resid_pre, resid_prb]), nanstd([resid_pre, resid_prb])./sqrt(sum(~isnan([resid_pre, resid_prb]))), 'k.-');

[h_res,p_res] = ttest(resid_pre - resid_prb);

%%
T_dev_fail = table; 
temp_dev = nan(50*num_locs*n_subs*(grp_days(i_grp) - 2)*2, 1);
temp_sub = nan(size(temp_dev));
temp_len = nan(size(temp_dev));
temp_prb = nan(size(temp_dev));
temp_day = nan(size(temp_dev));
i_ind = 1;
for i_sub = 1:n_subs
    for i_day = 3:grp_days(i_grp)
        temp_res = reshape(zu_rev(:, :, i_sub, i_day), 50, num_locs);
        temp_res = temp_res(sum(isnan(temp_res),2) < size(temp_res,2), :);
        temp_res = reshape(temp_res', numel(temp_res), 1);
        temp_dev(i_ind - 1 + (1:length(temp_res))) = temp_res;
        
        temp_sub(i_ind - 1 + (1:length(temp_res))) = i_sub;
        
        temp_len(i_ind - 1 + (1:length(temp_res))) = repmat(linspace(.25, .45, num_locs)', numel(temp_res)/num_locs, 1);
        
        temp_prb(i_ind - 1 + (1:length(temp_res))) = 0;
        
        temp_day(i_ind - 1 + (1:length(temp_res))) = i_day;
        i_ind = i_ind + length(temp_res);
        
        temp_res = reshape(zu_rev_probe(:, :, i_sub, i_day), 50, num_locs);
        temp_res = temp_res(sum(isnan(temp_res),2) < size(temp_res,2), :);
        temp_res = reshape(temp_res', numel(temp_res), 1);
        temp_dev(i_ind - 1 + (1:length(temp_res))) = temp_res;
        
        temp_sub(i_ind - 1 + (1:length(temp_res))) = i_sub;
        
        temp_len(i_ind - 1 + (1:length(temp_res))) = repmat(linspace(.25, .45, num_locs)', numel(temp_res)/num_locs, 1);
        
        if i_day == day_probe(i_grp)
            temp_prb(i_ind - 1 + (1:length(temp_res))) = 1;
        else
            temp_prb(i_ind - 1 + (1:length(temp_res))) = 0;
        end
        
        temp_day(i_ind - 1 + (1:length(temp_res))) = i_day;
        i_ind = i_ind + length(temp_res);
    end
end
T_dev_fail.len = double(temp_len);
T_dev_fail.sub = categorical(temp_sub);
T_dev_fail.dev = double(temp_dev);
T_dev_fail.prb = categorical(temp_prb);
T_dev_fail.succ = categorical(zeros(size(temp_len)));
T_dev_fail.day = categorical(temp_day);

%%
% score_by_trial_pre = reshape(nanmean(zu_rev, 2), 50, 20, grp_days(i_grp));
% score_by_trial_prb = reshape(nanmean(zu_rev_probe, 2), 50, 20, grp_days(i_grp));