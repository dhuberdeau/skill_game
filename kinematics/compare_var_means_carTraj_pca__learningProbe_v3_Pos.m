% for figure 6
%% load data

% load('traj_dir_15_90.mat')

% load('traj_dir_90.mat')
% load('tilt_dir_90.mat');
% traj_dir = tilt_dir;
load('traj_pos_25_90.mat');
load('exact_track_dist_full_v3.mat')
% load('exact_track_dist_full_strictExclusions_v2.mat') %changed to above
% 9/27/17

track_len_conv_fct = 226.8773;
for i_grp = 1:4
    traj_dir{i_grp} = traj_pos{i_grp}./track_len_conv_fct;
end

%% cancel subjects if needed:
sub_grp = [4];
sub_num = [12];

for i_sub = 1:length(sub_num)
    traj_dir{sub_grp(i_sub)}(:,:,sub_num(i_sub)) = nan;
end
%% setup parameters, variables, and constants

early_wind = {1:50, 201:250, 401:450, 601:650, 801:850, 1001:1050, 1201:1250, 1401:1450, 1601:1650, 1801:1850};
late_wind = {151:200, 351:400, 551:600, 751:800, 951:1000, 1151:1200, 1351:1400, 1551:1600, 1751:1800, 1951:2000};
pre_wind = {26:75, 226:275, 426:475, 626:675, 826:875, 1026:1075, 1226:1275, 1426:1475, 1626:1675, 1826:1875};
probe_wind = {76:125, 276:325, 476:525, 676:725, 876:925, 1076:1125, 1276:1325, 1476:1525, 1676:1725, 1876:1925};

other_wind = {...
    setdiff(26:200, [pre_wind{1} probe_wind{1}]),...1
    setdiff(226:400, [pre_wind{2} probe_wind{2}]),...2
    setdiff(426:600, [pre_wind{3} probe_wind{3}]),...3
    setdiff(626:800, [pre_wind{4} probe_wind{4}]),...4
    setdiff(826:1000, [pre_wind{5} probe_wind{5}]),...5
    setdiff(1026:1200, [pre_wind{6} probe_wind{6}]),...6
    setdiff(1226:1400, [pre_wind{7} probe_wind{7}]),...7
    setdiff(1426:1600, [pre_wind{8} probe_wind{8}]),...8
    setdiff(1626:1800, [pre_wind{9} probe_wind{9}]),...9
    setdiff(1826:2000, [pre_wind{10} probe_wind{10}])};%10

traj_var = {nan(84, 10, 4), nan(84, 10, 4), nan(84, 10, 4)}; %{all, succ, fail}(sub, day, pre/probe window)
traj_mdist = {nan(84, 10, 4), nan(84, 10, 4) nan(84, 10, 4)};
traj_kurt = {nan(84, 10, 4), nan(84, 10, 4), nan(84, 10, 4)}; %kurtosis
dist_trv = {nan(84, 10, 4), nan(84, 10, 4), nan(84, 10, 4)}; %distanct travelled, for purposes of correlation
traj_n = {nan(84, 10, 4), nan(84, 10, 4) nan(84, 10, 4)};
traj_covNorm = {nan(84, 10, 4), nan(84, 10, 4) nan(84, 10, 4)}; %covariance matrix norm
traj_mahal = {nan(84, 10, 4), nan(84, 10, 4), nan(84, 10, 4)};
% traj_basis = cell(1,10);

traj_pcData_pre = {nan(84, 10, 50), nan(84, 10, 50), nan(84, 10, 50)};
traj_pcData_prb = {nan(84, 10, 50), nan(84, 10, 50), nan(84, 10, 50)};

grp_s_skip = [0, 21, 42, 63];

loop_discard_TH = 200;
min_samp_n = 8;

% analyze only first part of signal:

% n_samp = 90;

k_samples = 1:60;
% k_samples = 21:80;

% which PC to use?
K_PC = 1;

%% collect distance travelled measures

day_group_inclusion = {1:4, 2:4, 2:4, 2:4, 2:4, 4, 4, 4, 4, 4};
for i_day = 1:10
    for i_grp = day_group_inclusion{i_day}
        for i_sub = 1:size(exact_track_dist_fall{i_grp},2)
            all_ind = early_wind{i_day};
            succ_ind = all_ind(~isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            fail_ind = all_ind(isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            dist_trv{1}(grp_s_skip(i_grp) + i_sub, i_day, 1) = nanmean(exact_track_dist_fall{i_grp}(all_ind, i_sub),1);
%             dist_trv{2}(grp_s_skip(i_grp) + i_sub, i_day, 1) = nanmean(exact_track_dist_fall{i_grp}(succ_ind, i_sub),1);
%             dist_trv{3}(grp_s_skip(i_grp) + i_sub, i_day, 1) = nanmean(exact_track_dist_fall{i_grp}(fail_ind, i_sub),1);
            
            all_ind = pre_wind{i_day};
            succ_ind = all_ind(~isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            fail_ind = all_ind(isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            dist_trv{1}(grp_s_skip(i_grp) + i_sub, i_day, 2) = nanmean(exact_track_dist_fall{i_grp}(all_ind, i_sub),1);
%             dist_trv{2}(grp_s_skip(i_grp) + i_sub, i_day, 2) = nanmean(exact_track_dist_fall{i_grp}(succ_ind, i_sub),1);
%             dist_trv{3}(grp_s_skip(i_grp) + i_sub, i_day, 2) = nanmean(exact_track_dist_fall{i_grp}(fail_ind, i_sub),1);
            
            all_ind = probe_wind{i_day};
            succ_ind = all_ind(~isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            fail_ind = all_ind(isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            dist_trv{1}(grp_s_skip(i_grp) + i_sub, i_day, 3) = nanmean(exact_track_dist_fall{i_grp}(all_ind, i_sub),1);
%             dist_trv{2}(grp_s_skip(i_grp) + i_sub, i_day, 3) = nanmean(exact_track_dist_fall{i_grp}(succ_ind, i_sub),1);
%             dist_trv{3}(grp_s_skip(i_grp) + i_sub, i_day, 3) = nanmean(exact_track_dist_fall{i_grp}(fail_ind, i_sub),1);
            
            all_ind = late_wind{i_day};
            succ_ind = all_ind(~isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            fail_ind = all_ind(isnan(traj_dir{i_grp}(end, all_ind, i_sub)));
            dist_trv{1}(grp_s_skip(i_grp) + i_sub, i_day, 4) = nanmean(exact_track_dist_fall{i_grp}(all_ind, i_sub),1);
%             dist_trv{2}(grp_s_skip(i_grp) + i_sub, i_day, 4) = nanmean(exact_track_dist_fall{i_grp}(succ_ind, i_sub),1);
%             dist_trv{3}(grp_s_skip(i_grp) + i_sub, i_day, 4) = nanmean(exact_track_dist_fall{i_grp}(fail_ind, i_sub),1);
        end
    end
end

%% Analyze car trajectory signals: PCA using uber-successful trajectories as basis for each subject
% uber-successful trajectories are ones passing 0.6
success_TH = .5;
signal_len = size(traj_dir{4},1)/2;
group_day_inclusion = {1, 1:5, 1:5, 1:10};
for i_grp = 1:length(traj_dir)    
    for i_sub = 1:size(traj_dir{i_grp},3)
        % Collect succ. traj.s OUTSIDE analysis regions for basis:
        traj_basis = nan(2*length(k_samples), 200*length(group_day_inclusion{i_grp})); %maximum size this matrix might be.. prune later
        k_basis = 0;
        for i_day = group_day_inclusion{i_grp}
            temp_base_a = traj_dir{i_grp}(k_samples, other_wind{i_day}, i_sub);
            temp_base_b = traj_dir{i_grp}(signal_len + k_samples, other_wind{i_day}, i_sub);
            temp_base = [temp_base_a; temp_base_b];
            
            temp_dist = exact_track_dist_fall{i_grp}(other_wind{i_day}, i_sub)';
            valid_inds = temp_dist >= success_TH; %only uber-successes.
            
            temp_base_2 = temp_base(:, valid_inds); 
            
            traj_basis(:, k_basis+(1:size(temp_base_2,2))) = temp_base_2;
            k_basis = k_basis + size(temp_base_2,2);
        end
        traj_basis = traj_basis(:, sum(isnan(traj_basis),1) < 1); %prune
        traj_basis_mean = nanmean(traj_basis, 2);
        traj_basis_ = traj_basis - repmat(traj_basis_mean, 1, size(traj_basis,2));
        if size(traj_basis_,2) > min_samp_n+1
            [traj_u, traj_s ,traj_v] = svd(traj_basis_, 'econ');        
            [jnk, jnk2, jnk3, jnk4, traj_explained] = pca(traj_basis_, 'Economy', 'on');

            for inclusion_criteria_state = 1
                for i_day = group_day_inclusion{i_grp}
                    temppre_1_a = traj_dir{i_grp}(k_samples, pre_wind{i_day}, i_sub);
                    temppre_1_b = traj_dir{i_grp}(signal_len + k_samples, pre_wind{i_day}, i_sub);
                    all_temppre_1_a = traj_dir{i_grp}(1:signal_len, pre_wind{i_day}, i_sub);
                    all_temppre_1_b = traj_dir{i_grp}(signal_len + (1:signal_len), pre_wind{i_day}, i_sub);
                    
                    switch inclusion_criteria_state
                        case 1
                            valid_inds = sum(isnan(temppre_1_a),1) < 1; % analyse ALL trajectories (mixure distribution)
                        case 2
                            valid_inds = sum(isnan(temppre_1_a),1) < 1 & ~isnan(all_temppre_1_a(end,:)); % analyse SUCCESSES only
                        case 3
                            valid_inds = sum(isnan(temppre_1_a),1) < 1 & isnan(all_temppre_1_a(end,:)); % analyse FAILURES only
                        otherwise
                            error('inclusion criteria mismatch');
                    end
                    temppre = [temppre_1_a(:, valid_inds); temppre_1_b(:, valid_inds)];
%                     temppre_3 = rad2deg(unwrap(deg2rad(temppre_2)));
%                     temppre = temppre_3(:, sum(abs(temppre_3) > loop_discard_TH, 1) < 1); %remove loops

                    tempp_1_a = traj_dir{i_grp}(k_samples, probe_wind{i_day}, i_sub);
                    tempp_1_b = traj_dir{i_grp}(signal_len + k_samples, probe_wind{i_day}, i_sub);
                    all_tempp_1_a = traj_dir{i_grp}(1:signal_len, probe_wind{i_day}, i_sub);
                    all_tempp_1_b = traj_dir{i_grp}(signal_len + (1:signal_len), probe_wind{i_day}, i_sub);
                    
                    switch inclusion_criteria_state
                        case 1
                            valid_inds = sum(isnan(tempp_1_a),1) < 1; % analyse ALL trajectories (mixure distribution)
                        case 2
                            valid_inds = sum(isnan(tempp_1_a),1) < 1 & ~isnan(all_tempp_1_a(end,:)); % analyse SUCCESSES only
                        case 3
                            valid_inds = sum(isnan(tempp_1_a),1) < 1 & isnan(all_tempp_1_a(end,:)); % analyse FAILURES only
                        otherwise
                            error('inclusion criteria mismatch');
                    end
                    tempp = [tempp_1_a(:, valid_inds); tempp_1_b(:, valid_inds)];
%                     tempp_3 = rad2deg(unwrap(deg2rad(tempp_2)));
%                     tempp = [tempp_1_a(:, sum(abs(tempp_3) > loop_discard_TH, 1) < 1); %remove loops

                    cond1 = size(temppre, 2) >= min_samp_n; cond2 = size(tempp, 2) >= min_samp_n; 

                    if cond1 && cond2
                        temppre_ = temppre - repmat(traj_basis_mean, 1, size(temppre,2));
                        tempp_ = tempp - repmat(traj_basis_mean, 1, size(tempp,2));

                        pca_pre = traj_u(:, K_PC)'*temppre_;
                        pca_all_pre = traj_u(:, :)'*temppre_;
%                         traj_var{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = std(pca_pre);
%                         traj_mdist{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = mean(pca_pre);
                        traj_var{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = std(pca_pre);
                        traj_mdist{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = mean(pca_pre);
                        traj_kurt{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = kurtosis(pca_pre);
                        traj_n{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = size(pca_pre,2);
                        traj_covNorm{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = norm(cov(pca_all_pre'));
                        for i_tr = 1:length(pca_pre) 
                            traj_pcData_pre{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, i_tr) = pca_pre(i_tr);
                        end
                        temp_mahal = nan(1, size(pca_all_pre, 2));
                        for i_tr = 1:length(temp_mahal)
                            temp_mahal(i_tr) = norm(pca_all_pre(:, i_tr));
                        end
                        traj_mahal{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 2) = nanmean(temp_mahal);
                        
                        pca_p = traj_u(:, K_PC)'*tempp_;
                        pca_all_p = traj_u(:, :)'*tempp_;
%                         traj_var{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = trace(sqrt(cov(pca_p)))/length(K_PC);
%                         traj_mdist{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = mean(pca_p(1,:));
                        traj_var{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = std(pca_p);
                        traj_mdist{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = mean(pca_p);
                        traj_kurt{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = kurtosis(pca_p(1,:));
                        traj_n{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = size(pca_p,2);
                        traj_covNorm{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = norm(cov(pca_all_p'));
                        for i_tr = 1:length(pca_p) 
                            traj_pcData_prb{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, i_tr) = pca_p(i_tr);
                        end
                        temp_mahal = nan(1, size(pca_all_p, 2));
                        for i_tr = 1:length(temp_mahal)
                            temp_mahal(i_tr) = norm(pca_all_p(:, i_tr));
                        end
                        traj_mahal{inclusion_criteria_state}(grp_s_skip(i_grp)+i_sub, i_day, 3) = nanmean(temp_mahal);
                    end
                end 
            end
        else
            warning(['No samples at criteria for this subject: ', num2str(i_sub)]);
        end
    end
end

%%

INCLUSION_STATE = 1; %1=mixture, 2=SUCC, 3=FAIL
state_clr = {'k', 'g', 'r'};

grp_sub_inds = {1:21, 22:42, 43:63, 64:84};
grp_all_day = {1, 1:5, 1:5, 1:10};
grp_train_day = {[], [1, 2, 4, 5], 1:4, 1:9};
grp_probe_day = {1, 3, 5, 10};
grp_1_inds = {1, 1:2:10, 1:2:10, 1:2:20};
grp_2_inds = {[], [2 4 8 10], [2 4 6 8], [2 4 6 8 10 12 14 16 18]};
grp_p_inds = [1, 6, 10, 20];
grp_x_train = {[.375 .625],...
    [.375 .625, 1.375 1.625, 2.375 2.625, 3.375 3.625, 4.375 4.625],...
    [.375 .625, 1.375 1.625, 2.375 2.625, 3.375 3.625, 4.375 4.625], ...
    [0.375 0.625, 1.375 1.625, 2.375 2.625, 3.375 3.625, 4.375 4.625, ...
    5.375 5.625, 6.375 6.625, 7.375 7.625, 8.375 8.625, 9.375 9.625]};
grp_x_probe = {.625, 2.625, 4.625, 9.625};
grp_name = {1,3,5,10};

min_ = nan; max_ = nan;
for i_grp = 2:4
    var_1 = traj_var{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 2);
    var_2 = traj_var{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 3);
    
    v_max = max(nanmean(var_1,1) + nanstd(var_1,0,1)./sqrt(sum(~isnan(var_1),1)));
    v_min = min(nanmean(var_1,1) - nanstd(var_1,0,1)./sqrt(sum(~isnan(var_1),1)));
    
    if isnan(min_)
        min_ = v_min;
    else
        if v_min < min_
            min_ = v_min;
        end
    end
    
    if isnan(max_)
        max_ = v_max;
    else
        if v_max > max_
            max_ = v_max;
        end
    end
end


% if ~exist('h1')
    h1 = figure;
% else
%     figure(h1);
% end
for i_grp = 2:4 %skip group 1, not enough data
    subplot(2,2,i_grp); hold on;
    var_1 = traj_var{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 2);
    var_2 = traj_var{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_train_day{i_grp}, 3);
    var_p = traj_var{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_probe_day{i_grp}, 3);
    
    y_train = nan(length(grp_sub_inds{i_grp}), length(grp_all_day{i_grp})*2);
    y_train(:, grp_1_inds{i_grp}) = var_1;
    y_train(:, grp_2_inds{i_grp}) = var_2;
    y_probe(:, 1) =  var_p;
    
    errorbar(grp_x_train{i_grp}, nanmean(y_train,1), nanstd(y_train,0,1)./sqrt(sum(~isnan(y_train),1)), [state_clr{INCLUSION_STATE}, '.-']);
    errorbar(grp_x_probe{i_grp}, nanmean(y_probe,1), nanstd(y_probe,0,1)./sqrt(sum(~isnan(y_probe),1)), [state_clr{INCLUSION_STATE}, 'o']);
    axis([0, 10, min_ - min_*.01, max_ + max_*.01]);
% axis([0, 10, 15, 25]);
    title(['Mixture Var, D', num2str(grp_name{i_grp})]);
end


min_ = 0; max_ = 0;
for i_grp = 2:4
    var_1 = traj_mdist{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 2);
    var_2 = traj_mdist{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 3);
    
    v_max = max(nanmean(var_1,1) + nanstd(var_1,0,1)./sqrt(sum(~isnan(var_1),1)));
    v_min = min(nanmean(var_1,1) - nanstd(var_1,0,1)./sqrt(sum(~isnan(var_1),1)));
    
    if isnan(min_)
        min_ = v_min;
    else
        if v_min < min_
            min_ = v_min;
        end
    end
    
    if isnan(max_)
        max_ = v_max;
    else
        if v_max > max_
            max_ = v_max;
        end
    end
end

% if ~exist('h2')
    h2 = figure;
% else
%     figure(h2);
% end
for i_grp = 2:4 %skip group 1, not enough data
    subplot(2,2,i_grp); hold on;
    d1_mean = nanmean(traj_mdist{1}(grp_sub_inds{i_grp}, 1, 2), 1);
    var_1 = traj_mdist{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_all_day{i_grp}, 2)-d1_mean;
    var_2 = traj_mdist{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_train_day{i_grp}, 3)-d1_mean;
    var_p = traj_mdist{INCLUSION_STATE}(grp_sub_inds{i_grp}, grp_probe_day{i_grp}, 3)-d1_mean;
    
    y_train = nan(length(grp_sub_inds{i_grp}), length(grp_all_day{i_grp})*2);
    y_train(:, grp_1_inds{i_grp}) = var_1;
    y_train(:, grp_2_inds{i_grp}) = var_2;
    y_probe(:, 1) =  var_p;
        
    errorbar(grp_x_train{i_grp}, nanmean(y_train,1), nanstd(y_train,0,1)./sqrt(sum(~isnan(y_train),1)), [state_clr{INCLUSION_STATE}, '.-']);
    errorbar(grp_x_probe{i_grp}, nanmean(y_probe,1), nanstd(y_probe,0,1)./sqrt(sum(~isnan(y_probe),1)), [state_clr{INCLUSION_STATE}, 'o']);
    axis([0, 10, min_ - min_*.01, max_ + max_*.01]);
% axis([0, 10, -15, 10]);
title(['Mixture Mean, D', num2str(grp_name{i_grp})]);
end
    
%% plot pooled curves 

sub_inds = 1:84;
day_probe_subs = {1:21, [], 22:42, [], 43:63, [], [], [], [], 64:84};

h3 = figure; hold on;

for INCLUSION_STATE = 1
    train_data = nan(84, 20);
    train_dist = nan(84, 20);
    probe_data = nan(21, 20);
    probe_dist = nan(21, 20);
    i_wind = 1;
    for i_day = 1:10
        daten_1 = traj_var{INCLUSION_STATE}(:, i_day, 2);
        daten_2 = traj_var{INCLUSION_STATE}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3);
        train_data(1:length(daten_1), i_wind) = daten_1;
        train_data(1:length(daten_2), i_wind + 1) = daten_2;

        daten_d1 = dist_trv{INCLUSION_STATE}(:, i_day, 2);
        daten_d2 = dist_trv{INCLUSION_STATE}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3);
        train_dist(1:length(daten_d1), i_wind) = daten_d1;
        train_dist(1:length(daten_d2), i_wind + 1) = daten_d2;

        daten_3 = traj_var{INCLUSION_STATE}(day_probe_subs{i_day}, i_day, 3);
        probe_data(1:length(daten_3), i_wind + 1) = daten_3;

        daten_d3 = dist_trv{INCLUSION_STATE}(day_probe_subs{i_day}, i_day, 3);
        probe_dist(1:length(daten_d3), i_wind + 1) = daten_d3;

        i_wind = i_wind + 2;
    end

    errorbar(grp_x_train{4}, nanmean(train_data, 1), nanstd(train_data)./sqrt(sum(~isnan(train_data),1)), [state_clr{INCLUSION_STATE}, '.-']);
    errorbar(grp_x_train{4}, nanmean(probe_data, 1), nanstd(probe_data)./sqrt(sum(~isnan(probe_data),1)), [state_clr{INCLUSION_STATE}, 'o']);
    % plot(repmat(grp_x_train{4}, size(train_data,1), 1), train_data, '.', 'Color', [.5 .5 .5]);
    % plot(repmat(grp_x_train{4}, size(probe_data,1), 1), probe_data, 'r.');
end

% h3_d = figure; hold on;
% plot(train_data(:), train_dist(:), '.');
% plot(probe_data(:), probe_dist(:), 'r.');
% title('Variance vs. distance travelled');
%%
h4 = figure; hold on;

for INCLUSION_STATE = 1
    train_data = nan(84, 20);
    train_dist = nan(84, 20);
    probe_data = nan(21, 20);
    probe_dist = nan(21, 20);
    
    train_frac = nan(84, 20);
    probe_frac = nan(21, 20);
    
    
    i_wind = 1;
    for i_day = 1:10
        daten_1 = traj_mdist{INCLUSION_STATE}(:, i_day, 2);
        daten_2 = traj_mdist{INCLUSION_STATE}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3);
        train_data(1:length(daten_1), i_wind) = daten_1;
        train_data(1:length(daten_2), i_wind + 1) = daten_2;

        daten_d1 = dist_trv{INCLUSION_STATE}(:, i_day, 2);
        daten_d2 = dist_trv{INCLUSION_STATE}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3);
        train_dist(1:length(daten_d1), i_wind) = daten_d1;
        train_dist(1:length(daten_d2), i_wind + 1) = daten_d2;
        
        daten_f1 = traj_n{2}(:, i_day, 2)./traj_n{1}(:, i_day, 2);
        daten_f2 = traj_n{2}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3)./traj_n{1}(setdiff(sub_inds, day_probe_subs{i_day}), i_day, 3);
        train_frac(1:length(daten_f1), i_wind) = daten_f1;
        train_frac(1:length(daten_f2), i_wind + 1) = daten_f2;

        daten_3 = traj_mdist{INCLUSION_STATE}(day_probe_subs{i_day}, i_day, 3);
        probe_data(1:length(daten_3), i_wind + 1) = daten_3;

        daten_d3 = dist_trv{INCLUSION_STATE}(day_probe_subs{i_day}, i_day, 3);
        probe_dist(1:length(daten_d3), i_wind + 1) = daten_d3;
        
        daten_f3 = traj_n{2}(day_probe_subs{i_day}, i_day, 3)./traj_n{1}(day_probe_subs{i_day}, i_day, 3);
        probe_frac(1:length(daten_d3), i_wind + 1) = daten_f3;

        i_wind = i_wind + 2;
    end

    errorbar(grp_x_train{4}, nanmean(train_data, 1), nanstd(train_data)./sqrt(sum(~isnan(train_data),1)), [state_clr{INCLUSION_STATE}, '.-']);
    errorbar(grp_x_train{4}, nanmean(probe_data, 1), nanstd(probe_data)./sqrt(sum(~isnan(probe_data),1)), [state_clr{INCLUSION_STATE}, 'o']);
    % plot(repmat(grp_x_train{4}, size(train_data,1), 1), train_data, '.', 'Color', [.5 .5 .5]);
    % plot(repmat(grp_x_train{4}, size(probe_data,1), 1), probe_data, 'r.');
end

% h_frac = figure; hold on;
% errorbar(grp_x_train{4}, nanmean(train_frac, 1), nanstd(train_frac)./sqrt(sum(~isnan(train_frac),1)), [state_clr{1}, '.-']);
% errorbar(grp_x_train{4}, nanmean(probe_frac, 1), nanstd(probe_frac)./sqrt(sum(~isnan(probe_frac),1)), [state_clr{1}, 'o']);

% h4_d = figure; hold on;
% plot(train_data(:), train_dist(:), '.');
% plot(probe_data(:), probe_dist(:), 'r.');
% title('Mean vs. distance travelled');

%% statistical analysis of probe variance:
model_single_var = cell(1,4);

% grp 2
i_grp = 2;
subjects_labs = reshape(repmat((1:21)', 1, 5), 21*5, 1);
train_labs = zeros(21*5,1);
probe_labs = [zeros(21*2, 1); ones(21, 1); zeros(21*2, 1)];
pre_wind_wind = reshape(repmat([.375 1.375, 2.375, 3.375 4.375], 21, 1), 21*5, 1);
prb_wind_wind = reshape(repmat([.625 1.625, 2.625, 3.625 4.625], 21, 1), 21*5, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_var{INCLUSION_STATE}(22:42, 1:5, 2), 21*5, 1);...
    reshape(traj_var{INCLUSION_STATE}(22:42, 1:5, 3), 21*5, 1)];

T_var_2 = table;
T_var_2.subject = [subjects_labs; subjects_labs];
T_var_2.probe = [train_labs; probe_labs];
T_var_2.window = window_labs;
T_var_2.response = datenmatrix;
try
lm_var_2 = fitlme(T_var_2, 'response ~ window + probe + (1|subject)');
model_single_var{i_grp} = lm_var_2;
catch 
end


% grp 3
i_grp = 3;
subjects_labs = reshape(repmat((1:21)', 1, 5), 21*5, 1);
train_labs = zeros(21*5,1);
probe_labs = [zeros(21*4, 1); ones(21, 1)];
pre_wind_wind = reshape(repmat([.375 1.375, 2.375, 3.375 4.375], 21, 1), 21*5, 1);
prb_wind_wind = reshape(repmat([.625 1.625, 2.625, 3.625 4.625], 21, 1), 21*5, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_var{INCLUSION_STATE}(43:63, 1:5, 2), 21*5, 1);...
    reshape(traj_var{INCLUSION_STATE}(43:63, 1:5, 3), 21*5, 1)];

T_var_3 = table;
T_var_3.subject = [subjects_labs; subjects_labs];
T_var_3.probe = [train_labs; probe_labs];
T_var_3.window = window_labs;
T_var_3.response = datenmatrix;
try
lm_var_3 = fitlme(T_var_3, 'response ~ window + probe + (1|subject)');
model_single_var{i_grp} = lm_var_3;
catch
end

% grp 4
i_grp = 4;
subjects_labs = reshape(repmat((1:21)', 1, 10), 21*10, 1);
train_labs = zeros(21*10,1);
probe_labs = [zeros(21*9, 1); ones(21, 1)];
pre_wind_wind = reshape(repmat([[.375 1.375, 2.375, 3.375 4.375], 5+[.375 1.375, 2.375, 3.375 4.375]], 21, 1), 21*10, 1);
prb_wind_wind = reshape(repmat([[.625 1.625, 2.625, 3.625 4.625], 5+[.625 1.625, 2.625, 3.625 4.625]], 21, 1), 21*10, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_var{INCLUSION_STATE}(64:84, 1:10, 2), 21*10, 1);...
    reshape(traj_var{INCLUSION_STATE}(64:84, 1:10, 3), 21*10, 1)];

T_var_4 = table;
T_var_4.subject = [subjects_labs; subjects_labs];
T_var_4.probe = [train_labs; probe_labs];
T_var_4.window = window_labs;
T_var_4.response = datenmatrix;
try
lm_var_4 = fitlme(T_var_4, 'response ~ window + probe + (1|subject)');
model_single_var{i_grp} = lm_var_4;
catch
end
%% statistical analysis of probe means:
model_single_mean = cell(1,4);

% % grp 1
% i_grp = 1;
% subjects_labs = reshape(repmat((1:21)', 1, 5), 21*5, 1);
% train_labs = zeros(21*5,1);
% probe_labs = [zeros(21*2, 1); ones(21, 1); zeros(21*2, 1)];
% pre_wind_wind = reshape(repmat([.375 1.375, 2.375, 3.375 4.375], 21, 1), 21*5, 1);
% prb_wind_wind = reshape(repmat([.625 1.625, 2.625, 3.625 4.625], 21, 1), 21*5, 1);
% window_labs = [pre_wind_wind; prb_wind_wind];
% datenmatrix = [reshape(traj_mdist{INCLUSION_STATE}(22:42, 1:5, 2), 21*5, 1);...
%     reshape(traj_mdist{INCLUSION_STATE}(22:42, 1:5, 3), 21*5, 1)];
% 
% T_var_2 = table;
% T_var_2.subject = [subjects_labs; subjects_labs];
% T_var_2.probe = [train_labs; probe_labs];
% T_var_2.window = window_labs;
% T_var_2.response = datenmatrix;
% lm_var_2 = fitlme(T_var_2, 'response ~ window + probe + (1|subject)');
% model_single_mean{i_grp} = lm_var_2;

% grp 2
i_grp = 2;
subjects_labs = reshape(repmat((1:21)', 1, 5), 21*5, 1);
train_labs = zeros(21*5,1);
probe_labs = [zeros(21*2, 1); ones(21, 1); zeros(21*2, 1)];
pre_wind_wind = reshape(repmat([.375 1.375, 2.375, 3.375 4.375], 21, 1), 21*5, 1);
prb_wind_wind = reshape(repmat([.625 1.625, 2.625, 3.625 4.625], 21, 1), 21*5, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_mdist{INCLUSION_STATE}(22:42, 1:5, 2), 21*5, 1);...
    reshape(traj_mdist{INCLUSION_STATE}(22:42, 1:5, 3), 21*5, 1)];

T_var_2 = table;
T_var_2.subject = [subjects_labs; subjects_labs];
T_var_2.probe = [train_labs; probe_labs];
T_var_2.window = window_labs;
T_var_2.response = datenmatrix;
lm_var_2 = fitlme(T_var_2, 'response ~ window + probe + (1|subject)');
model_single_mean{i_grp} = lm_var_2;

% grp 3
i_grp = 3;
subjects_labs = reshape(repmat((1:21)', 1, 5), 21*5, 1);
train_labs = zeros(21*5,1);
probe_labs = [zeros(21*4, 1); ones(21, 1)];
pre_wind_wind = reshape(repmat([.375 1.375, 2.375, 3.375 4.375], 21, 1), 21*5, 1);
prb_wind_wind = reshape(repmat([.625 1.625, 2.625, 3.625 4.625], 21, 1), 21*5, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_mdist{INCLUSION_STATE}(43:63, 1:5, 2), 21*5, 1);...
    reshape(traj_mdist{INCLUSION_STATE}(43:63, 1:5, 3), 21*5, 1)];

T_var_3 = table;
T_var_3.subject = [subjects_labs; subjects_labs];
T_var_3.probe = [train_labs; probe_labs];
T_var_3.window = window_labs;
T_var_3.response = datenmatrix;
lm_var_3 = fitlme(T_var_3, 'response ~ window + probe + (1|subject)');
model_single_mean{i_grp} = lm_var_3;

% grp 4
i_grp = 4;
subjects_labs = reshape(repmat((1:21)', 1, 10), 21*10, 1);
train_labs = zeros(21*10,1);
probe_labs = [zeros(21*9, 1); ones(21, 1)];
pre_wind_wind = reshape(repmat([[.375 1.375, 2.375, 3.375 4.375], 5+[.375 1.375, 2.375, 3.375 4.375]], 21, 1), 21*10, 1);
prb_wind_wind = reshape(repmat([[.625 1.625, 2.625, 3.625 4.625], 5+[.625 1.625, 2.625, 3.625 4.625]], 21, 1), 21*10, 1);
window_labs = [pre_wind_wind; prb_wind_wind];
datenmatrix = [reshape(traj_mdist{INCLUSION_STATE}(64:84, 1:10, 2), 21*10, 1);...
    reshape(traj_mdist{INCLUSION_STATE}(64:84, 1:10, 3), 21*10, 1)];

T_var_4 = table;
T_var_4.subject = [subjects_labs; subjects_labs];
T_var_4.probe = [train_labs; probe_labs];
T_var_4.window = window_labs;
T_var_4.response = datenmatrix;
lm_var_4 = fitlme(T_var_4, 'response ~ window + probe + (1|subject)');
model_single_mean{i_grp} = lm_var_4;

%%

df3 = traj_var{1}(22:42, 3,  3) - traj_var{1}(22:42, 3, 2);
df5 = traj_var{1}(43:63, 5,  3) - traj_var{1}(43:63, 5, 2);
df10 = traj_var{1}(64:84, 10,  3) - traj_var{1}(64:84, 10, 2);

figure; 
subplot(1,2,1); hold on;
bar([1 2 3], [nanmean(df3), nanmean(df5), nanmean(df10)], 'k');
errorbar([1 2 3], [nanmean(df3), nanmean(df5), nanmean(df10)],...
    [nanstd(df3), nanstd(df5), nanstd(df10)]./sqrt(sum(~isnan([df3, df5, df10]),1)), 'k.');
axis([0 4 -.1 .1])
[h_var, p_var, junk, d_var] = ttest([df3, df5, df10]);
p_var_new = pval_adjust(p_var, 'holm')';
[~, p_var_all, ~, d_var_all] = ttest([df3; df5; df10]);

df3 = traj_mdist{1}(22:42, 3,  3) - traj_mdist{1}(22:42, 3, 2);
df5 = traj_mdist{1}(43:63, 5,  3) - traj_mdist{1}(43:63, 5, 2);
df10 = traj_mdist{1}(64:84, 10,  3) - traj_mdist{1}(64:84, 10, 2);

subplot(1,2,2); hold on;
bar([1 2 3], [nanmean(df3), nanmean(df5), nanmean(df10)], 'k');
errorbar([1 2 3], [nanmean(df3), nanmean(df5), nanmean(df10)],...
    [nanstd(df3), nanstd(df5), nanstd(df10)]./sqrt(sum(~isnan([df3, df5, df10]),1)), 'k.');
axis([0 4 -.01 .01])
[h_m, p_m, junk, d_m] = ttest([df3, df5, df10]);
p_m_new = pval_adjust(p_m, 'holm')';
[~, p_m_all, ~, d_m_all] = ttest([df3; df5; df10]);