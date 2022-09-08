library(latex2exp)

# Fit shifted logarithm and polynomial models to Aubert et al. telomere length 
#  data and recover division estimates

set.seed(123)
df <- read.csv("~/PathToData");

data <- df[,c("age", "Granulo")];
data <- na.omit(data);
x <- data[["age"]];
y <- data[["Granulo"]];

############################ poly-fit

poly_model <- lm(y ~ poly(x, 3, raw=T))
summary(poly_model)

poly_AIC = AIC(poly_model)

coef <- unname(poly_model$coefficients)

d_30 <- coef /-0.03
d_65 <- coef / -0.065
d_100 <- coef / -0.1

polyDiv <- function(t, params) {
  c <- params[2]
  b <- params[3]
  a <- params[4]
  predict <- a*t**3 + b*t**2 + c*t
}
plot(x,polyDiv(x,d_30), type="l", col="red", lwd=3, xlab="Age (years)", ylab="Divisions")
lines(x, polyDiv(x,d_65), col="black", lwd=5)
lines(x, polyDiv(x,d_100), col="blue", lwd=3)
grid(col="#5D6D7E")
legend(0,130, c("30 bp","65 bp", "100 bp"), lwd=c(3,5,3), col=c("red","black", "blue"), title=TeX(r'($\Delta L$)'), bg = "white")


########################## log-fit


log_model <- nls(y ~ a + b*log(c*x+1), start=c(a=10,b=-0.3, c=0.02), lower=c(5,-10,0), upper=c(15,0, 1), algorith="port" )

summary(log_model)

log_coef <- unname(coef(log_model))
log_AIC = AIC(log_model)

logDiv <- function(t, params) {
  a2 <- params[1]
  a3 <- params[2]
  predict <- a2*log(a3*t +1)
  
}

log_d_30 <- c(log_coef[2]/-0.03, log_coef[3])
log_d_65 <- c(log_coef[2]/-0.065, log_coef[3])
log_d_100 <- c(log_coef[2]/-0.1, log_coef[3])

plot(x,logDiv(x,log_d_30), type="l", col="red", lwd=3, xlab="Age (years)", ylab="Divisions")
lines(x, logDiv(x,log_d_65), col="black", lwd=5)
lines(x, logDiv(x,log_d_100), col="blue", lwd=3)

grid(col="#5D6D7E")
legend(0,120, c("30 bp","65 bp", "100 bp"), lwd=c(3,5,3), col=c("red","black", "blue"), title=TeX(r'($\Delta L$)'), bg = "white")