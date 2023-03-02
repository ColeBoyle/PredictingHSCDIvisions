library(segmented)

set.seed(123)
df <- read.csv("AndreuGranData.csv", na.strings="-");
data <- df[,c("Age","Gran")];
data <- na.omit(data);
x <- data[["Age"]];
y <- data[["Gran"]];


lf <- lm(y ~ x);
phase_1_AIC = AIC(lf)
phase_1_BIC = BIC(lf)
summary(lf)


################################ 2-phase

phase_2 <- segmented(lf, npsi=1)
summary(phase_2)

phase_2_AIC = AIC(phase_2)
phase_2_BIC = BIC(phase_2)

slope_2 <- slope(phase_3)

################################ 3-phase

phase_3 <- segmented(lf, npsi=2)
summary(phase_3)

slope_3 <- slope(phase_3)

phase_3_AIC = AIC(phase_3)
phase_3_BIC = BIC(phase_3)
################################ 4-phase

#phase_4 <- segmented.lm(lf, npsi=)
#summary(phase_4)
#
#phase_4_AIC = AIC(phase_4)
#phase_4_BIC = BIC(phase_4)
#slope(phase_4)
