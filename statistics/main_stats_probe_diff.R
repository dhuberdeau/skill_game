# Title     : TODO
# Objective : Test for changes in the difference in mean and 
#             variability (standard deviation) among trajectories
#             across groups, across principle component, and from 
#             baseline (intercept - which is possibly the most
#             important one to test to see if mean or var. 
#             is significantly different during probes.)
# Created by: david
# Created on: 9/4/18

# Load libraries 
library(pwr)

# Set working directory
# setwd("/Users/david/OneDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/public_analysis/statistics")
setwd("/Users/david/Desktop/iPad_analysis_temp/public_analysis/statistics")

# define function:
probe_diff_fcn <- function(mData_1){
  # Load libraries 
  
  mX_1 <- data.frame(y = mData_1$V1, group = factor(mData_1$V2), subj = mData_1$V4)
  mod_1 <- lm(y ~ group, data = mX_1)
  
  return(mod_1)
}

########
## Analyze distance travelled measures
# Load mean data (created from Matlab)
data_dist <- read.csv("../data/dist_travelled_probe_diff", header=FALSE)
mod_dist <- probe_diff_fcn(data_dist)
stats_dist <- summary(mod_dist)
p_adj_dist <- p.adjust(stats_dist$coefficients[,4], method = "holm")

########
## Analyze kinematics measures (mean and variability of PCs)
# Load mean data (created from Matlab)
data_mean <- read.csv("../data/mean_probe_difference_pc1", header=FALSE)
mod_mean <- probe_diff_fcn(data_mean)
stats_mean <- summary(mod_mean)
p_adj_mean <- p.adjust(stats_mean$coefficients[,4], method = "holm")

data_amean <- data_mean
data_amean$V1 <- abs(data_mean$V1)
mod_amean <- probe_diff_fcn(data_amean)
stats_amean <- summary(mod_amean)
p_adj_amean <- p.adjust(stats_amean$coefficients[,4], method = "holm")

# Load variance data (created from Matlab):
data_var <- read.csv("../data/var_probe_difference_pc1", header=FALSE)
mod_var <- probe_diff_fcn(data_var)
stats_var <- summary(mod_var)
p_adj_var <- p.adjust(stats_var$coefficients[,4], method = "holm")

########
## Analyze policy deviation measures
# Load mean data (created from Matlab)
sdata_policy <- read.csv("../data/policy_successes_probe_diff_grp", header=FALSE)
sdata_policy_ex <- subset(sdata_policy, V2 > 1, select = c(V1, V2, V3, V4))
smod_policy <- probe_diff_fcn(sdata_policy_ex)
sstats_policy <- summary(smod_policy)
sp_adj_policy <- p.adjust(sstats_policy$coefficients[,4], method = "holm")

fdata_policy <- read.csv("../data/policy_failures_probe_diff_grp", header=FALSE)
fmod_policy <- probe_diff_fcn(fdata_policy)
fstats_policy <- summary(fmod_policy)
fp_adj_policy <- p.adjust(fstats_policy$coefficients[,4], method = "holm")

cdata_policy <- read.csv("../data/policy_union_probe_diff_grp", header=FALSE)
cmod_policy <- probe_diff_fcn(cdata_policy)
cstats_policy <- summary(cmod_policy)
cp_adj_policy <- p.adjust(cstats_policy$coefficients[,4], method = "holm")


#######
## Analyze movement variability measures for difference between pre-probe and probe 
## lumping D3, D5, and D10 groups together in order to gain better power:

var_all <- data_var$V1[data_var$V2 > 1]
var_t_results <- t.test(var_all, mu = 0, alternative = "two.sided")
t_results <- var_t_results

########
# determine power for one-sample t-test (which is what we have used):
n_ = length(var_all)
es_ = .005 #t_results$estimate
se_ = sqrt(var(var_all)/length(var_all))
sd_ = se_*sqrt(n_)
d_ <- es_/sd_
pw_var = pwr.t.test(d = d_, sig.level = 0.05, n = n_, alternative = "two.sided", 
                  type = "one.sample")

#########
# analyze power of difference computation individually, using the approximate values below:
n_ = 20
es_ = mean(c(sstats_policy$coefficients[1], sstats_policy$coefficients[2], sstats_policy$coefficients[3]))
se_ = mean(c(sstats_policy$coefficients[4], sstats_policy$coefficients[5], sstats_policy$coefficients[6]))
sd_ = se_*sqrt(n_)
d_ <- es_/sd_
pw_ = pwr.t.test(d = d_, sig.level = 0.05, n = n_, alternative = "two.sided", 
                 type = "one.sample")

#######
## Analyze policy deviation measures for difference between pre-probe and probe 
## lumping D3, D5, and D10 groups together in order to gain better power:

s_policy_all <- sdata_policy$V1[sdata_policy$V2 > 1]
s_t_results <- t.test(s_policy_all, mu = 0, alternative = "two.sided")
t_results <- s_t_results

########
# determine power for one-sample t-test (which is what we have used):
n_ = length(s_policy_all)
es_ = .1 #t_results$estimate
se_ = sqrt(var(s_policy_all)/length(s_policy_all))
sd_ = se_*sqrt(n_)
d_ <- es_/sd_
pw_s = pwr.t.test(d = d_, sig.level = 0.05, n = n_, alternative = "two.sided", 
           type = "one.sample")

#######
## Analyze policy deviation measures for difference between pre-probe and probe 
## lumping D3, D5, and D10 groups together in order to gain better power:

f_policy_all <- fdata_policy$V1[fdata_policy$V2 > 1]
f_t_results <- t.test(f_policy_all, mu = 0, alternative = "two.sided")
t_results <- f_t_results

########
# determine power for one-sample t-test (which is what we have used):
n_ = length(f_policy_all)
es_ = .1 #t_results$estimate
se_ = sqrt(var(f_policy_all)/length(f_policy_all))
sd_ = se_*sqrt(n_)
d_ <- es_/sd_
pw_f = pwr.t.test(d = d_, sig.level = 0.05, n = n_, alternative = "two.sided", 
                 type = "one.sample")

# d = effect size / standard deviation