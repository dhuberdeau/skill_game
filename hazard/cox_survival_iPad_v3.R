library(survival)
library(coxme)
library(effects)


setwd("~/OneDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/complete_analysis")
#setwd("~/OneDrive/Documents/JHU/BLAM_lab/R_dir/iPad_skill")
#setwd("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/R_dir/iPad_skill") #Office PC

#E1_cat <- read.csv("event_category", header = FALSE, na.strings = 'NaN')
#E1_rate <- read.csv("hazard_response_night_early", header = FALSE, na.strings = 'NaN')
#E1_rate <- read.csv("event_times", header = FALSE, na.strings = 'NaN')
E1_cat <- read.csv("event_category_window_25_pre_prb2", header = FALSE, na.strings = 'NaN')
E1_rate <- read.csv("event_times_window_25_pre_prb2", header = FALSE, na.strings = 'NaN')

dat <- data.frame(subject = (E1_cat[,1]),
                  group = factor(E1_cat[,2]),
                  day = factor(E1_cat[,3]),
                  probe = factor(E1_cat[,4]),
                  probe2 = factor(E1_cat[,5]),
                  prePrb = factor(E1_cat[,6]),
                  window = factor(E1_cat[,7]),
                  events = E1_rate[,1])

#exclude D1
dat2 <- data.frame(subject = (E1_cat[4201:84200,1]),
                  group = factor(E1_cat[4201:84200,2]),
                  day = factor(E1_cat[4201:84200,3]),
                  probe = factor(E1_cat[4201:84200,4]),
                  probe2 = factor(E1_cat[4201:84200,5]),
                  prePrb = factor(E1_cat[4201:84200,6]),
                  window = factor(E1_cat[4201:84200,7]),
                  events = E1_rate[4201:84200,1])

#exclude D1, D10
dat3 <- data.frame(subject = (E1_cat[4201:44200,1]),
                   group = factor(E1_cat[4201:44200,2]),
                   day = factor(E1_cat[4201:44200,3]),
                   probe = factor(E1_cat[4201:44200,4]),
                   probe2 = factor(E1_cat[4201:44200,5]),
                   prePrb = factor(E1_cat[4201:44200,6]),
                   window = factor(E1_cat[4201:44200,7]),
                   events = E1_rate[4201:44200,1])

#exclude D1, D5
dat4 <- data.frame(subject = (E1_cat[c(4201:24200, 44201:84200),1]),
                   group = factor(E1_cat[c(4201:24200, 44201:84200),2]),
                   day = factor(E1_cat[c(4201:24200, 44201:84200),3]),
                   probe = factor(E1_cat[c(4201:24200, 44201:84200),4]),
                   probe2 = factor(E1_cat[c(4201:24200, 44201:84200),5]),
                   prePrb = factor(E1_cat[c(4201:24200, 44201:84200),6]),
                   window = factor(E1_cat[c(4201:24200, 44201:84200),7]),
                   events = E1_rate[c(4201:24200, 44201:84200),1])

#exclude D1, D3
dat5 <- data.frame(subject = (E1_cat[24201:84200,1]),
                   group = factor(E1_cat[24201:84200,2]),
                   day = factor(E1_cat[24201:84200,3]),
                   probe = factor(E1_cat[24201:84200,4]),
                   probe2 = factor(E1_cat[24201:84200,5]),
                   prePrb = factor(E1_cat[24201:84200,6]),
                   window = factor(E1_cat[24201:84200,7]),
                   events = E1_rate[24201:84200,1])

#Create the simplest test data set: each group individually

# Fit a stratified model for group 1
# m1_probe <- coxme(Surv(events) ~ probe + (1|subject), subset = group == 1 & 
#                     (window == 3 | window == 4 | window == 5), data = dat)
# m1_ <- coxme(Surv(events) ~ (1|subject), subset = group == 1 & (window == 3 | window == 4 | window == 5), data = dat)
# m1_probe <- coxme(Surv(events) ~ probe2 + (1|subject), subset = group == 1 & (day == 1), data = dat)
# m1_ <- coxme(Surv(events) ~ (1|subject), subset = group == 1 & day == 1, data = dat)
m1_probe <- coxme(Surv(events) ~ probe2 + (1|subject), 
                  subset = group == 1 & (probe2 == 1 | prePrb == 1), 
                  data = dat)
m1_ <- coxme(Surv(events) ~ (1|subject), 
             subset = group == 1 & (probe2 == 1 | prePrb == 1), 
             data = dat)
comp_probe1 <- anova(m1_, m1_probe)

# Fit a stratified model for group 3
# m2_probe <- coxme(Surv(events) ~ probe + (1|subject), subset = group == 2 & 
#                     (window == 19 | window == 20 | window == 21), data = dat)
# m2_ <- coxme(Surv(events) ~ (1|subject), subset = group == 2 & 
#                (window == 19 | window == 20 | window == 21), data = dat)
m2_probe <- coxme(Surv(events) ~ probe2 + (1|subject), 
                  subset = group == 2 & (probe2 == 1 | prePrb == 1),
                  data = dat)
m2_ <- coxme(Surv(events) ~ (1|subject), 
             subset = group == 2 & (probe2 == 1 | prePrb == 1),
             data = dat)
comp_probe2 <- anova(m2_, m2_probe)

# Fit a stratified model for group 5
# m3_probe <- coxme(Surv(events) ~ probe + (1|subject), subset = group == 3 & 
#                     (window == 35 | window == 36 | window == 37), data = dat)
# m3_ <- coxme(Surv(events) ~ (1|subject), subset = group == 3 & 
#                (window == 35 | window == 36 | window == 37), data = dat)
m3_probe <- coxme(Surv(events) ~ probe2 + (1|subject), 
                  subset = group == 3 & (probe2 == 1 | prePrb == 1), 
                  data = dat)
m3_ <- coxme(Surv(events) ~ (1|subject), 
             subset = group == 3 & (probe2 == 1 | prePrb == 1), 
             data = dat)
comp_probe3 <- anova(m3_, m3_probe)

# Fit a stratified model for group 10
# m4_probe <- coxme(Surv(events) ~ probe + (1|subject), subset = group == 4 & 
#                     (window == 75 | window == 76 | window == 77), data = dat)
# m4_ <- coxme(Surv(events) ~ (1|subject), subset = group == 4 & 
#                (window == 75 | window == 76 | window == 77), data = dat)
m4_probe <- coxme(Surv(events) ~ probe2 + (1|subject), 
                  subset = group == 4 & (probe2 == 1 | prePrb == 1), 
                  data = dat)
m4_ <- coxme(Surv(events) ~ (1|subject), 
             subset = group == 4 & (probe2 == 1 | prePrb == 1), 
             data = dat)
comp_probe4 <- anova(m4_, m4_probe)

# plot barchart with error bars of standard error
grp_coefs <- c(m1_probe$coefficients[1], m2_probe$coefficients[1], m3_probe$coefficients[1], m4_probe$coefficients[1])
grp_se <- c(0.044, 0.042, 0.037, 0.040)
barCenters <- barplot(grp_coefs, ylim = c(-0.3, 0.6))
arrows(barCenters,  grp_coefs - grp_se, barCenters, grp_coefs + grp_se, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

# # Test day, probe, and interaction across all groups
mall_x <- coxme(Surv(events) ~ probe2 + group + probe2:group + (1|subject), 
                subset = (probe2 == 1 | prePrb == 1) & (group == 1 & day == 1) | (group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10),
                data = dat)
mall_d_p <- coxme(Surv(events) ~ probe2 + group + (1|subject), 
                  subset = (probe2 == 1 | prePrb == 1) & (group == 1 & day == 1) | (group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10),
                  data = dat)
mall_p <- coxme(Surv(events) ~ probe2 + (1|subject), 
                subset = (probe2 == 1 | prePrb == 1) & (group == 1 & day == 1) | (group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10),
                data = dat)
mall_d <- coxme(Surv(events) ~ group + (1|subject), 
                subset = (probe2 == 1 | prePrb == 1) & (group == 1 & day == 1) | (group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10),
                data = dat)

cx <- anova(mall_x, mall_d_p)
cp <- anova(mall_d_p, mall_p)
cd <- anova(mall_d_p, mall_d)
# 
# # Test day, probe, and interaction across only groups 3, 5, & 10
# # Fit day model and test against no factors

mall_x3 <- coxme(Surv(events) ~ probe2 + group + probe2:group + (1|subject), 
                 subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10)),
                 data = dat2)
mall_d_p3 <- coxme(Surv(events) ~ probe + group + (1|subject), 
                   subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10)),
                  data = dat2)
mall_p3 <- coxme(Surv(events) ~ probe + (1|subject), 
                 subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10)),
                data = dat2)
mall_d3 <- coxme(Surv(events) ~ group + (1|subject), 
                 subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5) | (group == 4 & day == 10)),
                data = dat2)

cx3 <- anova(mall_x3, mall_d_p3)
cp3 <- anova(mall_d_p3, mall_p3)
cd3 <- anova(mall_d_p3, mall_d3)

# Test day, probe, and interaction between 3 & 5, 3 & 10, and 5 & 10
mall_x3_5 <- coxme(Surv(events) ~ group + probe2 + group:probe2 + (1|subject), 
                   subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5)),
                   data = dat3)
mall_p_d3_5 <- coxme(Surv(events) ~ group + probe2 + (1|subject), 
                     subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) |(group == 3 & day == 5)),
                     data = dat3)
mall_p3_5 <- coxme(Surv(events) ~ probe2 + (1|subject), 
                   subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5)),
                   data = dat3)
mall_d3_5 <- coxme(Surv(events) ~ group + (1|subject), 
                   subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 3 & day == 5)),
                   data = dat3)

cx3_5 <- anova(mall_x3_5, mall_p_d3_5)
cp3_5 <- anova(mall_p_d3_5, mall_p3_5)
cd3_5 <- anova(mall_p_d3_5, mall_d3_5)

mall_x3_10 <- coxme(Surv(events) ~ probe2*group + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 4 & day == 10)), data = dat4)
mall_d_p3_10 <- coxme(Surv(events) ~ group + probe2 + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 4 & day == 10)), data = dat4)
mall_p3_10 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 4 & day == 10)), data = dat4)
mall_d3_10 <- coxme(Surv(events) ~ group + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 2 & day == 3) | (group == 4 & day == 10)), data = dat4)

cx3_10 <- anova(mall_x3_10, mall_d_p3_10)
cp3_10 <- anova(mall_d_p3_10, mall_p3_10)
cd3_10 <- anova(mall_d_p3_10, mall_d3_10)

mall_x5_10 <- coxme(Surv(events) ~ probe2*group + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 3 & day == 5) | (group == 4 & day == 10)), data = dat5)
mall_d_p5_10 <- coxme(Surv(events) ~ group + probe2 + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 3 & day == 5) | (group == 4 & day == 10)), data = dat5)
mall_p5_10 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 3 & day == 5) | (group == 4 & day == 10)), data = dat5)
mall_d5_10 <- coxme(Surv(events) ~ group + (1|subject), subset = (probe2 == 1 | prePrb == 1) & ((group == 3 & day == 5) | (group == 4 & day == 10)), data = dat5)

cx5_10 <- anova(mall_x5_10, mall_d_p5_10)
cp5_10 <- anova(mall_d_p5_10, mall_p5_10)
cd5_10 <- anova(mall_d_p5_10, mall_d5_10)
# 
# 
# ## plot stuff
# s_day <- survfit(Surv(events) ~ strata(day), data = dat)
# #plot(s_day)
# 
# clr_scale <- gray.colors(10, start = .7, end = 0, gamma = 2.2, alpha = NULL);
# plot(s_day, col = clr_scale)
# 
# s_probe <- survfit(Surv(events) ~ strata(probe), data = dat)
# plot(s_probe)
# # test and plot pre-probe vs. probe:
s_prb1 <- survfit(Surv(events) ~ 1, subset = group == 1 & day == 1 & prePrb == 1, data = dat)
s_prb1_p <- survfit(Surv(events) ~ 1, subset = group == 1 & day == 1 & probe2 == 1, data = dat)
plot(s_prb1, col = "darkblue")
lines(s_prb1_p, col = "brown1")
m_prb1 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = day == 1 & (probe2 == 1 | prePrb == 1), data = dat)
# 
s_prb2 <- survfit(Surv(events) ~ 1, subset = group == 2 & day == 3 & prePrb == 1, data = dat)
s_prb2_p <- survfit(Surv(events) ~ 1, subset = group == 2 & day == 3 & probe2 == 1, data = dat)
plot(s_prb2, col = "darkblue")
lines(s_prb2_p, col = "brown1")
m_prb2 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = day == 3 & (probe2 == 1 | prePrb == 1), data = dat)
# 
s_prb3 <- survfit(Surv(events) ~ strata(probe2), subset = group == 3 & day == 5 & prePrb == 1, data = dat)
s_prb3_p <- survfit(Surv(events) ~ strata(probe2), subset = group == 3 & day == 5 & probe2 == 1, data = dat)
plot(s_prb3, col = "darkblue")
lines(s_prb3_p, col = "brown1")
m_prb3 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = day == 5 & (probe2 == 1 | prePrb == 1), data = dat)
# 
s_prb4 <- survfit(Surv(events) ~ strata(probe2), subset = group == 4 & day == 10 & prePrb == 1, data = dat)
s_prb4_p <- survfit(Surv(events) ~ strata(probe2), subset = group == 4 & day == 10 & probe2 == 1, data = dat)
plot(s_prb4, col = "darkblue")
lines(s_prb4_p, col = "brown1")
m_prb4 <- coxme(Surv(events) ~ probe2 + (1|subject), subset = day == 10 & (probe2 == 1 | prePrb == 1), data = dat)

n1 <- 1
dt <- 1
while(dt > 0) {dt <- .5 - s_prb1$time[n1]; n1 <- n1 + 1;}

n2 <- 1
dt <- 1
while(dt > 0) {dt <- .5 - s_prb1$time[n2]; n2 <- n2 + 1;}

n3 <- 1
dt <- 1
while(dt > 0) {dt <- .5 - s_prb1$time[n3]; n3 <- n3 + 1;}

n4 <- 1
dt <- 1
while(dt > 0) {dt <- .5 - s_prb1$time[n4]; n4 <- n4 + 1;}

grp_crossSec <- c(s_prb1$surv[n1] - s_prb1_p$surv[n1], 
                  s_prb2$surv[n2] - s_prb2_p$surv[n2], 
                  s_prb3$surv[n3] - s_prb3_p$surv[n3], 
                  s_prb4$surv[n4] - s_prb4_p$surv[n4])

grp_cs_se <- c(s_prb1$std.err[n1], s_prb2$std.err[n2], s_prb3$std.err[n3], s_prb4$std.err[n4])

barCenters <- barplot(grp_crossSec, ylim = c(-0.3, 0.6))
arrows(barCenters,  grp_crossSec - grp_cs_se, barCenters, grp_crossSec + grp_cs_se, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


plot(s_prb1$time, s_prb1$surv, type = "l", col = "orange")
lines(s_prb1$time, s_prb1$surv + s_prb1$std.err, col = "orange")
lines(s_prb1$time, s_prb1$surv - s_prb1$std.err, col = "orange")
lines(s_prb1_p$time, s_prb1_p$surv, type = "l", col = "black")
lines(s_prb1_p$time, s_prb1_p$surv + s_prb1_p$std.err, col = "black")
lines(s_prb1_p$time, s_prb1_p$surv - s_prb1_p$std.err, col = "black")

plot(s_prb2$time, s_prb2$surv, type = "l", col = "red")
lines(s_prb2$time, s_prb2$surv + s_prb2$std.err, col = "red")
lines(s_prb2$time, s_prb2$surv - s_prb2$std.err, col = "red")
lines(s_prb2_p$time, s_prb2_p$surv, col = "black")
lines(s_prb2_p$time, s_prb2_p$surv + s_prb2_p$std.err, col = "black")
lines(s_prb2_p$time, s_prb2_p$surv - s_prb2_p$std.err, col = "black")

plot(s_prb3$time, s_prb3$surv, type = "l", col = "blue")
lines(s_prb3$time, s_prb3$surv - s_prb3$std.err, col = "blue")
lines(s_prb3$time, s_prb3$surv + s_prb3$std.err, col = "blue")
lines(s_prb3_p$time, s_prb3_p$surv, col = "black")
lines(s_prb3_p$time, s_prb3_p$surv - s_prb3_p$std.err, col = "black")
lines(s_prb3_p$time, s_prb3_p$surv + s_prb3_p$std.err, col = "black")

plot(s_prb4$time, s_prb4$surv, type = "l", col = "green")
lines(s_prb4$time, s_prb4$surv + s_prb4$std.err, col = "green")
lines(s_prb4$time, s_prb4$surv - s_prb4$std.err, col = "green")
lines(s_prb4_p$time, s_prb4_p$surv, col = "black")
lines(s_prb4_p$time, s_prb4_p$surv + s_prb4_p$std.err, col = "black")
lines(s_prb4_p$time, s_prb4_p$surv - s_prb4_p$std.err, col = "black")