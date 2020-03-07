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
probe_wind_fcn <- function(mData_1){
  # Load libraries 
  # library(lme4)
  
  mX_1 <- data.frame(y = mData_1$V1, wind = mData_1$V2, group = factor(mData_1$V3), subj = mData_1$V4)
  mod_1 <- lm(y ~ 1 + wind, data = mX_1)
  
  return(mod_1)
}

compute_standardized_betas <- function(mod){
  # compute the estimate's effect size relative to 
  library("QuantPsyc")
  coef_lmbeta <- lm.beta(mod)
  return(coef_lmbeta)
}

########
## Analyze distance travelled measures
# Load mean data (created from Matlab)
data_dist <- read.csv("../data/dist_travelled_window", header=FALSE)
dist_1 <- subset(data_dist, V3 == 1, select = c(V1, V2, V3, V4))
dist_2 <- subset(data_dist, V3 == 2, select = c(V1, V2, V3, V4))
dist_3 <- subset(data_dist, V3 == 3, select = c(V1, V2, V3, V4))
dist_4 <- subset(data_dist, V3 == 4, select = c(V1, V2, V3, V4))
mod_dist <- probe_wind_fcn(data_dist)
mod_dist_1 <- probe_wind_fcn(dist_1)
mod_dist_2 <- probe_wind_fcn(dist_2)
mod_dist_3 <- probe_wind_fcn(dist_3)
mod_dist_4 <- probe_wind_fcn(dist_4)
stats_dist <- summary(mod_dist)
p_adj_dist <- p.adjust(stats_dist$coefficients[,4], method = "holm")

########
## Analyze kinematics measures (mean and variability of PCs)
# Load mean data (created from Matlab)
data_mean <- read.csv("../data/mean_data_pc1", header=FALSE)
mean_1 <- subset(data_mean, V3 == 1, select = c(V1, V2, V3, V4))
#mean_1$V1 <- 2^(abs(mean_1$V1))
#mean_1$V2 <- mean_1$V2 - 100
mean_2 <- subset(data_mean, V3 == 2, select = c(V1, V2, V3, V4))
#mean_2$V1 <- 2^(abs(mean_2$V1))
#mean_2$V2 <- mean_2$V2 - 500
mean_3 <- subset(data_mean, V3 == 3, select = c(V1, V2, V3, V4))
#mean_3$V1 <- 2^(abs(mean_3$V1))
#mean_3$V2 <- mean_3$V2 - 500
mean_4 <- subset(data_mean, V3 == 4, select = c(V1, V2, V3, V4))
#mean_4$V1 <- 2^(abs(mean_4$V1))
#mean_4$V2 <- mean_4$V2 - 1000
mod_mean <- probe_wind_fcn(data_dist)

mod_mean_1 <- probe_wind_fcn(mean_1)
mbeta_1 <- compute_standardized_betas(mod_mean_1)
mod_mean_2 <- probe_wind_fcn(mean_2)
mbeta_2 <- compute_standardized_betas(mod_mean_2)
mod_mean_3 <- probe_wind_fcn(mean_3)
mbeta_3 <- compute_standardized_betas(mod_mean_3)
mod_mean_4 <- probe_wind_fcn(mean_4)
mbeta_4 <- compute_standardized_betas(mod_mean_4)

mod_mean <- probe_wind_fcn(data_mean)
stats_mean <- summary(mod_mean)
p_adj_mean <- p.adjust(stats_mean$coefficients[,4], method = "holm")

# Load mean absolute value data (created from Matlab)
data_mean <- read.csv("../data/mean_data_pc1", header=FALSE)
data_mean$V1 <- exp(abs(data_mean$V1))
amean_1 <- subset(data_mean, V3 == 1, select = c(V1, V2, V3, V4))
#mean_1$V1 <- 2^(abs(mean_1$V1))
#mean_1$V2 <- mean_1$V2 - 100
amean_2 <- subset(data_mean, V3 == 2, select = c(V1, V2, V3, V4))
#mean_2$V1 <- 2^(abs(mean_2$V1))
#mean_2$V2 <- mean_2$V2 - 500
amean_3 <- subset(data_mean, V3 == 3, select = c(V1, V2, V3, V4))
#mean_3$V1 <- 2^(abs(mean_3$V1))
#mean_3$V2 <- mean_3$V2 - 500
amean_4 <- subset(data_mean, V3 == 4, select = c(V1, V2, V3, V4))
#mean_4$V1 <- 2^(abs(mean_4$V1))
#mean_4$V2 <- mean_4$V2 - 1000
mod_mean <- probe_wind_fcn(data_dist)

mod_amean_1 <- probe_wind_fcn(amean_1)
ambeta_1 <- compute_standardized_betas(mod_amean_1)
mod_amean_2 <- probe_wind_fcn(amean_2)
ambeta_2 <- compute_standardized_betas(mod_amean_2)
mod_amean_3 <- probe_wind_fcn(amean_3)
ambeta_3 <- compute_standardized_betas(mod_amean_3)
mod_amean_4 <- probe_wind_fcn(amean_4)
ambeta_4 <- compute_standardized_betas(mod_amean_4)

mod_amean <- probe_wind_fcn(data_mean)
stats_amean <- summary(mod_amean)
p_adj_amean <- p.adjust(stats_amean$coefficients[,4], method = "holm")


# Load variance data (created from Matlab):
data_var <- read.csv("../data/variability_data_pc1", header=FALSE)
data_var$V1 <- log(data_var$V1) # log transform the responses
mod_var <- probe_wind_fcn(data_var)
var_1 <- subset(data_var, V3 == 1, select = c(V1, V2, V3, V4))
var_2 <- subset(data_var, V3 == 2, select = c(V1, V2, V3, V4))
var_3 <- subset(data_var, V3 == 3, select = c(V1, V2, V3, V4))
var_4 <- subset(data_var, V3 == 4, select = c(V1, V2, V3, V4))
mod_var <- probe_wind_fcn(data_var)

mod_var_1 <- probe_wind_fcn(var_1)
vbeta_1 <- compute_standardized_betas(mod_var_1)
mod_var_2 <- probe_wind_fcn(var_2)
vbeta_2 <- compute_standardized_betas(mod_var_2)
mod_var_3 <- probe_wind_fcn(var_3)
vbeta_3 <- compute_standardized_betas(mod_var_3)
mod_var_4 <- probe_wind_fcn(var_4)
vbeta_4 <- compute_standardized_betas(mod_var_4)

stats_var <- summary(mod_var)
p_adj_var <- p.adjust(stats_var$coefficients[,4], method = "holm")

########
## Analyze policy deviation measures
# Load mean data (created from Matlab)
data_policy <- read.csv("../data/policy_successes", header=FALSE)
smod_policy <- probe_wind_fcn(data_policy)
spolicy_1 <- subset(data_policy, V3 == 1, select = c(V1, V2, V3, V4))
spolicy_2 <- subset(data_policy, V3 == 2, select = c(V1, V2, V3, V4))
spolicy_3 <- subset(data_policy, V3 == 3, select = c(V1, V2, V3, V4))
spolicy_4 <- subset(data_policy, V3 == 4, select = c(V1, V2, V3, V4))
smod_policy <- probe_wind_fcn(data_policy)
smod_policy_1 <- probe_wind_fcn(spolicy_1)
smod_policy_2 <- probe_wind_fcn(spolicy_2)
smod_policy_3 <- probe_wind_fcn(spolicy_3)
smod_policy_4 <- probe_wind_fcn(spolicy_4)
sstats_policy <- summary(smod_policy)
sp_adj_policy <- p.adjust(sstats_policy$coefficients[,4], method = "holm")

data_policy <- read.csv("../data/policy_failures", header=FALSE)
smod_policy <- probe_wind_fcn(data_policy)
fpolicy_1 <- subset(data_policy, V3 == 1, select = c(V1, V2, V3, V4))
fpolicy_2 <- subset(data_policy, V3 == 2, select = c(V1, V2, V3, V4))
fpolicy_3 <- subset(data_policy, V3 == 3, select = c(V1, V2, V3, V4))
fpolicy_4 <- subset(data_policy, V3 == 4, select = c(V1, V2, V3, V4))
fmod_policy <- probe_wind_fcn(data_policy)
fmod_policy_1 <- probe_wind_fcn(fpolicy_1)
fmod_policy_2 <- probe_wind_fcn(fpolicy_2)
fmod_policy_3 <- probe_wind_fcn(fpolicy_3)
fmod_policy_4 <- probe_wind_fcn(fpolicy_4)
fstats_policy <- summary(fmod_policy)
fp_adj_policy <- p.adjust(fstats_policy$coefficients[,4], method = "holm")

data_policy <- read.csv("../data/policy_union", header=FALSE)
cmod_policy <- probe_wind_fcn(data_policy)
cpolicy_1 <- subset(data_policy, V3 == 1, select = c(V1, V2, V3, V4))
cpolicy_2 <- subset(data_policy, V3 == 2, select = c(V1, V2, V3, V4))
cpolicy_3 <- subset(data_policy, V3 == 3, select = c(V1, V2, V3, V4))
cpolicy_4 <- subset(data_policy, V3 == 4, select = c(V1, V2, V3, V4))
cmod_policy <- probe_wind_fcn(data_policy)
cmod_policy_1 <- probe_wind_fcn(cpolicy_1)
cmod_policy_2 <- probe_wind_fcn(cpolicy_2)
cmod_policy_3 <- probe_wind_fcn(cpolicy_3)
cmod_policy_4 <- probe_wind_fcn(cpolicy_4)
cstats_policy <- summary(cmod_policy)
cp_adj_policy <- p.adjust(cstats_policy$coefficients[,4], method = "holm")


# determine power for one-sample t-test (which is what we have used):
d <- 0.50/2.25
pwr.t.test(d = d, sig.level = 0.05, n = 20, alternative = "two.sided", 
           type = "one.sample")

# d = effect size / standard deviation


# determine power (for linear model):
u_ = length(fstats_policy$terms)
n_ = 81
r_ = fstats_policy$r.squared
sl_ = .05
f_pow = pwr.f2.test(u = u_, v = n_ - u_ - 1, f2 = r_/(1 - r_), sig.level = sl_)

u_ = length(sstats_policy$terms)
n_ = 81
r_ = sstats_policy$r.squared
sl_ = .05
s_pow = pwr.f2.test(u = u_, v = n_ - u_ - 1, f2 = r_/(1 - r_), sig.level = sl_)

u_ = length(fstats_policy$terms)
n_ = 81
r_ = .1
sl_ = .05
g_pow = pwr.f2.test(u = u_, v = n_ - u_ - 1, f2 = r_/(1 - r_), sig.level = sl_)

##      Multiple regression power calculation 
## 
##               u = 2 (degrees of freedom = number of independent variables)
##               v = 49.88971 (degrees of freedom = n - u - 1)
##              f2 = 0.4285714 (effect size (r^2/(1-r^2)))
##       sig.level = 0.001 (usually .05)
##           power = 0.8 (desired power)