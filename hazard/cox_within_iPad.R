library(survival)
library(coxme)

# setwd("~/OneDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/complete_analysis")
setwd("~/OneDrive/Documents/JHU/BLAM_lab/R_dir/iPad_skill")

E1_cat <- read.csv("event_category_winday_25", header = FALSE, na.strings = 'NaN')
E1_rate <- read.csv("event_times_winday_25", header = FALSE, na.strings = 'NaN')

dat <- data.frame(subject = factor(E1_cat[,1]),
                  day = E1_cat[,2],
                  window = factor(E1_cat[,3]),
                  events = E1_rate[,1])

#Create the simplest test data set 

# Fit a models
m1 <- coxme(Surv(events) ~ day*window + (1|subject), data = dat) 
m_x <- coxme(Surv(events) ~ day + window + (1|subject), data = dat)
m_d <- coxme(Surv(events) ~ day + (1|subject), data = dat)
m_w <- coxme(Surv(events) ~ window + (1|subject), data = dat)

m_comp_x<-anova(m1, m_x)
m_comp_d<-anova(m_x, m_d)
m_comp_w<-anova(m_x, m_w)


s1_w1 <- survfit(Surv(events) ~ strata(day), subset = window == 0, data = dat)
s1_w2 <- survfit(Surv(events) ~ strata(day), subset = window == 1, data = dat)
plot(s1_w1, col = "black")
lines(s1_w2, col = "red")

w0_temp <- summary(s1_w1, times = 0.5)
w1_temp <- summary(s1_w2, times = 0.5)
x_temp <- c(.3, 1.3, 2.3, 3.3, 4.3, 5.3, 6.3, 7.3, 8.3, 9.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
plot(x = x_temp, y = c(w0_temp$surv, w1_temp$surv), type = "p", asp = 5)
arrows(x_temp, c(w0_temp$surv - w0_temp$std.err, w1_temp$surv - w1_temp$std.err),
       x_temp, c(w0_temp$surv + w0_temp$std.err, w1_temp$surv + w1_temp$std.err),
       length = 0.05, angle = 90, code = 3)
for(i in 1:10){
  lines(x=c(x_temp[i], x_temp[i+10]), y=c(w0_temp$surv[i], w1_temp$surv[i]), asp = 5)
}


s1_wa <- survfit(Surv(events) ~ 1, subset = window == 0, data = dat)
s1_wb <- survfit(Surv(events) ~ 1, subset = window == 1, data = dat)
plot(s1_wa, col = "darkblue")
lines(s1_wb, col = "brown1")

