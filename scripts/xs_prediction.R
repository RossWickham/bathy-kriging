library(reshape2)
library(rgl)
library(tidyverse)

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
points(data$qout, fit$fitted.values, col="red", pch=3, cex=0.1)
summary(fit)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.839e-03  5.674e-04   3.242   0.0012 ** 
#   poly(qout, 2, raw = T)1 -7.312e-08  1.013e-08  -7.219 5.99e-13 ***
#   poly(qout, 2, raw = T)2  1.094e-11  3.329e-14 328.583  < 2e-16 ***
#   ---

data$pred <- coef(fit)[1] + data$qout*coef(fit)[2] + data$qout^2*coef(fit)[3]

points(data$qout, data$pred, col = "green")

# Predicting past data
load("c:/gitrepo/bathy-kriging/data/cwms_data_2020-2025.RData")

dts <- intersect(export_data$datetime[export_data$path == "IHR.Elev-Forebay.Inst.1Hour.0.CBT-REV"],
                 export_data$datetime[export_data$path == "LMN.Flow-Out.Ave.1Hour.1Hour.CBT-REV"])

export_data %>%
  split(export_data$path) %>%
  lapply(nrow)

cwmsData <- export_data %>%
  mutate(
    path = str_replace_all(path, "IHR.Elev-Forebay.Inst.1Hour.0.CBT-REV", "fb"),
    path = str_replace_all(path, "LMN.Flow-Out.Ave.1Hour.1Hour.CBT-REV" , "qout"),
    dateValue = as.numeric(datetime)
  ) %>%
  filter(dateValue %in% dts) %>%
  arrange(dateValue) %>%
  na.omit() %>%
  dcast(formula = "datetime ~ path", value.var = "value", fun.aggregate = max) %>%
  na.omit()


