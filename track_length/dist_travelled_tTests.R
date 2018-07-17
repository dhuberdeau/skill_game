library(lme4)
library(boot)
library(effects)

#E1_rate <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E1_rate", header = FALSE, na.strings = 'na')
#E1_cat <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E1_categorical", header = FALSE, na.strings = 'na')

#E1_rate <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/withinDay_data", na.strings = 'NaN')
#E1_cat <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/withinDay_design", na.strings = 'NaN')
#p_values <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/mdt_diff_p", header = FALSE, na.strings = 'NaN')
p_values <- read.csv("/Users/hal/Desktop/SkyDrive/Documents/JHU/BLAM_lab/Projects/skillLearning_antGame_A/Matlab/nouveau_analysis/analysis_scripts/distance_travel_v2/haz_diff_p", header = FALSE, na.strings = 'NaN')

#E1_rate <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_rate_xcat2", header = FALSE, na.strings = 'NaN')
#E1_cat <- read.csv("C:/Users/David/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_categorical_xcat2", header = FALSE, na.strings = 'NaN')


#E2_rate <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_rate", header = FALSE, colClasses = "numeric")
#E2_cat <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E2_categorical", header = FALSE, colClasses = "factor")

#E3_rate <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E3_rate", header = FALSE, colClasses = "numeric")
#E3_cat <- read.csv("C:/Users/JH/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_timedResponse/meta_analysis/E3_categorical", header = FALSE, colClasses = "factor")


p_adj <- p.adjust(p_values, method = p.adjust.methods[3], n = 4)

summary(p_adj)