library(rgl)


load("M:\\CRT\\XS_Stage_Prediction\\Analysis Results\\Year_Round\\R06\\Results_DailyMax_2017-09-29_01_32_25\\R06 Regression Results.RData")


plot3d(allData$Lag000_ICE_HARBOR_POOL,
       allData$Lag000_LOWER_MONUMENTAL_OUT,
       allData$SNAKE_REACH_06__24.1755- allData$Lag000_ICE_HARBOR_POOL)


data <- data.frame(fb = allData$Lag000_ICE_HARBOR_POOL,
                   qout = allData$Lag000_LOWER_MONUMENTAL_OUT,
                   wse = allData$SNAKE_REACH_06__24.1755)

data$wse_delta <- data$wse - data$fb

fit <- lm(wse_delta ~ poly(qout, 2, raw=T), data)

plot(data$qout, data$wse_delta, pch=3, cex=0.1)
points(data$qout, fit$fitted.values, col="red")
summary(fit)

data$pred <- coef(fit)[1] + data$qout*coef(fit)[2] + data$qout^2*coef(fit)[3]

points(data$qout, data$pred, col = "green")


