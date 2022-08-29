library(segmented)

# Fit 2-4 phase piecewise linear models to Aubert et al. telomere length data

set.seed(123)
df <- read.csv("~/Path/To/AubertData");

data <- df[,c("age", "Granulo")];
data <- na.omit(data);
x <- data[["age"]];
y <- data[["Granulo"]];

lf <- lm(y ~ x);
################################ 2-phase

phase_2 <- segmented.lm(lf, npsi=1)
summary(phase_2)

phase_2_AIC = AIC(phase_2)
slope(phase_2) 
################################ 3-phase

phase_3 <- segmented.lm(lf, npsi=2)
summary(phase_3)

phase_3_AIC = AIC(phase_3)
slope(phase_3)
################################ 4-phase

phase_4 <- segmented.lm(lf, npsi=3)
summary(phase_4)

phase_4_AIC = AIC(phase_4)
slope(phase_4)