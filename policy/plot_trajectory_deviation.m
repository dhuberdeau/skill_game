

% len_range = [.25 .7];
TRS_DAY = 200;
TRS_WIND = 50;
grp_days = [1 5 5 10];
n_subs = [21 20 20 20];
% num_locs = round(100*(len_range(2) - len_range(1)));
num_locs = ceil(100*(len_range(2) - len_range(1)));
pre_inds = {26:75, 426:475, 826:875, 1826:1875};
probe_inds = {76:125, 476:525, 876:925, 1876:1925};
pre_inds_day = repmat((26:75)', 1, 10) + repmat(0:200:1999, 50, 1);
probe_inds_day = repmat((101:125)', 1, 10) + repmat(0:200:1999, 25, 1);

zu_trial_pre_successes = cell(1,4);
for i_grp = 1:4
    zu_trial_pre_successes{i_grp} = ...
        nan(size(pre_inds_day,1), num_locs, n_subs(i_grp), grp_days(i_grp));
    for i_day = 1:grp_days(i_grp)
        temp_score_day_pre = score_s{i_grp}(pre_inds_day(:, i_day), :, :);
        zu_trial_pre_successes{i_grp}(:, :, :, i_day) = temp_score_day_pre;
    end
end

zu_trial_probe_successes = cell(1,4);
for i_grp = 1:4
    zu_trial_probe_successes{i_grp} = ...
        nan(size(probe_inds_day,1), num_locs, n_subs(i_grp), grp_days(i_grp));
    for i_day = 1:grp_days(i_grp)
        temp_score_day_pre = score_s{i_grp}(probe_inds_day(:, i_day), :, :);
        zu_trial_probe_successes{i_grp}(:, :, :, i_day) = temp_score_day_pre;
    end
end


REMOVE_SUB = [];

pre_wind = {26:75, 226:275, 426:475, 626:675, 826:875, 1026:1075, 1226:1275, 1426:1475, 1626:1675, 1826:1875};
probe_wind = {76:125, 276:325, 476:525, 676:725, 876:925, 1076:1125, 1276:1325, 1476:1525, 1676:1725, 1876:1925};

for i_grp = 1:4
    zu_trial = zu_trial_pre_successes{i_grp}(:, :, 1:n_subs(i_grp), 1:grp_days(i_grp));
    zu_trial_probe = zu_trial_probe_successes{i_grp}(:, :, 1:n_subs(i_grp), 1:grp_days(i_grp));

    zu_fwd = nan(size(zu_trial));
    zu_fwd_probe = nan(size(zu_trial_probe));
    zu_rev = nan(size(zu_trial));
    zu_rev_probe = nan(size(zu_trial_probe));
    loc_inds = 1:num_locs;
    for i_day = 1:grp_days(i_grp)
        for i_sub = setdiff(1:size(zu_trial,3), REMOVE_SUB)
            for i_tr = 1:size(zu_trial,1)
                temp = zu_trial(i_tr, :, i_sub, i_day);
                k_max = max(loc_inds(~isnan(temp)));
                temp2 = temp(1:k_max);
                zu_fwd(i_tr, 1:k_max, i_sub, i_day) = temp2;
                zu_rev(i_tr, (num_locs - k_max + 1):num_locs, i_sub, i_day) = temp2;
            end
            for i_tr = 1:size(zu_trial_probe,1)
                temp = zu_trial_probe(i_tr, :, i_sub, i_day);
                k_max = max(loc_inds(~isnan(temp)));
                temp2 = temp(1:k_max);
                zu_fwd_probe(i_tr, 1:k_max, i_sub, i_day) = temp2;
                zu_rev_probe(i_tr, (num_locs - k_max + 1):num_locs, i_sub, i_day) = temp2;
            end
        end
    end


    figure;
    for i_day = 1:grp_days(i_grp)
        subplot(2,5,i_day); hold on;
        dev_pre = nan(n_subs(i_grp), num_locs);
        dev_prb = nan(n_subs(i_grp), num_locs);
        for i_sub = setdiff(1:20, REMOVE_SUB)
            temp = zu_rev(:,:,i_sub, i_day);
            temp_pre = temp(sum(isnan(temp),2) < size(temp,2), :);

            temp = zu_rev_probe(:,:,i_sub, i_day);
            temp_prb = temp(sum(isnan(temp),2) < size(temp,2), :);

            dev_pre(i_sub, :) = nanmean(temp_pre,1);
            dev_prb(i_sub, :) = nanmean(temp_prb,1);
        end
        errorfield(1:num_locs, nanmean(dev_pre,1), nanstd(dev_pre)./sqrt(sum(~isnan(dev_pre(:,end)))), 'k');
        errorfield(1:num_locs, nanmean(dev_prb,1), nanstd(dev_prb)./sqrt(sum(~isnan(dev_prb(:,end)))), 'r');
        axis([0 num_locs 0 2])
        title('Reverse');
    end
end