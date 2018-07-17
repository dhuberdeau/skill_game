library(lme4)
library(boot)
library(effects)

#E1_rate <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E1_rate", header = FALSE, na.strings = 'na')
#E1_cat <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E1_categorical", header = FALSE, na.strings = 'na')

#E1_rate <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/withinDay_data", na.strings = 'NaN')
#E1_cat <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/withinDay_design", na.strings = 'NaN')
#E1_rate <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/acrossDay_data", na.strings = 'NaN')
#E1_cat <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/acrossDay_design", na.strings = 'NaN')

#E1_rate <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_rate_xcat2", header = FALSE, na.strings = 'NaN')
#E1_cat <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_categorical_xcat2", header = FALSE, na.strings = 'NaN')


#E2_rate <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_rate", header = FALSE, colClasses = "numeric")
#E2_cat <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_categorical", header = FALSE, colClasses = "factor")

E1_rate <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E3_rate", header = FALSE, colClasses = "numeric")
E1_cat <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E3_categorical", header = FALSE, colClasses = "factor")

dat <- data.frame(subject = E1_cat[,1],
                  Y1 = E1_rate[,1])

dat$subject <- factor(dat$subject)
dat$Y1 <- as.numeric(dat$Y1)

lm1 <- lmer(Y1 ~ (1|subject), data = dat, REML = FALSE)

summary(lm1)
#lm2 <- update(lm1, formula = . ~ . + type:cycle)
#lm3 <- lmer(Y1 ~ type + (1+cycle+type|subject), data = dat, REML = FALSE)
#lm4 <- lmer(Y1 ~ cycle + (1+cycle+type|subject), data = dat, REML = FALSE)

#anova(lm1, lm2) # LRT
#anova(lm1, lm3)
#anova(lm1, lm4)
#plot(allEffects(lm2))
# an1 <- aov(Y1 ~ cycle * type, data = dat)
# m1 <- lm(Y1 ~ cycle * type, data = dat)

#tmpfun <- function(.) {
#  fxf <- allEffects(.)$type$fit
#  c(fxf[4:6] - fxf[1:3])
#}
#b1 <- bootMer(lm2, tmpfun, nsim = 100, seed = 1, 
#.progress = 'txt')
#?bootMer