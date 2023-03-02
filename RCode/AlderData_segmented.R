library(segmented)
set.seed(123)

Alder_gran_data <- c()
    
Alder_age_data <- c()

data <- data.frame(Alder_age_data, Alder_gran_data)
data <- na.omit(data) 
data <- data[order(data$Alder_age_data),]
x <-  data[["Alder_age_data"]]
y <-  data[["Alder_gran_data"]]


lf <- lm(y ~ x);
summary(lf)
phase_1_AIC = AIC(lf)
phase_1_BIC = BIC(lf)
################################ 2-phase

phase_2 <- segmented(lf, npsi =1)
summary(phase_2)

slope_2 = slope(phase_2)
phase_2_AIC = AIC(phase_2)
phase_2_BIC = BIC(phase_2)

################################ 3-phase

phase_3 <- segmented(lf, npsi =2)
summary(phase_3)

################################ 4-phase

phase_4 <- segmented(lf, npsi=3)
summary(phase_4)

phase_4_AIC = AIC(phase_4)
phase_4_BIC = BIC(phase_4)
slope(phase_4)


################################ 6-phase

#phase_5 <- segmented(lf, npsi=4)
#summary(phase_5)
#
#phase_5_AIC = AIC(phase_5)
#phase_5_BIC = BIC(phase_5)
#slope(phase_5)

