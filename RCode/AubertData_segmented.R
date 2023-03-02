library(segmented)
# Fit 1-4 phase piecewise linear models to Aubert et al. telomere length data

set.seed(123)
df <- read.csv("PLOS-2012.csv")

data <- df[,c("age", "Granulo")];
data <- na.omit(data);
x <- data[["age"]];
y <- data[["Granulo"]];
lf <- lm(y ~ x);
summary(lf)
plot(y ~ x)

phase_1_AIC = AIC(lf)
phase_1_BIC = BIC(lf)

################################ 2-phase

phase_2 <- segmented(lf, npsi=1)
summary(phase_2)
slope(phase_2)
slope_2 <- slope(phase_2)

phase_2_AIC = AIC(phase_2)
phase_2_BIC = BIC(phase_2)
################################ 3-phase

phase_3 <- segmented(lf, npsi=2)
summary(phase_3)
slope_3 <- slope(phase_3)

phase_3_AIC = AIC(phase_3)
phase_3_BIC = BIC(phase_3)
################################ 4-phase

phase_4 <- segmented(lf, npsi=3)
summary(phase_4)
slope_4 <- slope(phase_4)
phase_4_AIC = AIC(phase_4)
phase_4_BIC = BIC(phase_4)


############################### 5-phase
#phase_5 <- segmented(lf, npsi=4)
#summary(phase_5)
#
#phase_5_AIC = AIC(phase_5)
#phase_5_BIC = BIC(phase_5)
#
#plot(x,y, col="grey")
#plot(phase_5, add=T)

############################### 6-phase
phase_6 <- segmented(lf, npsi=5)
summary(phase_6)

phase_6_AIC = AIC(phase_6)
phase_6_BIC = BIC(phase_6)


############################### 7-phase
phase_7 <- segmented(lf, npsi=6)
summary(phase_7)

phase_7_AIC = AIC(phase_7)
phase_7_BIC = BIC(phase_7)
