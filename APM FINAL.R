################################################################################
## APM - Assignment 2
################################################################################
# Load required packages -------------------------------------------------------
library(ggplot2)
library(lubridate)
library(forecast)
library(tseries)
library(fpp)
library(gridExtra)
library(astsa)
library(xts)
library(tidyverse)
library(rio)
library(reshape2)
library(ggcorrplot)
library(GGally)
library(TSstudio)
library(urca)

# Read the data ----------------------------------------------------------------
pollution <- read.csv("http://www.stats.gla.ac.uk/~tereza/rp/Air_pollution_assignment.csv")

# Convert Date column to POSIXct format
pollution$Date <- as.Date(data$Date, format = "%Y-%m-%d %H:%M:%S")

## EDA #########################################################################

# Check the structure of the data ----------------------------------------------
str(data)

ggplot(pollution, aes(Date, PM10)) +
  geom_line(color = "#41ab5d") +
  scale_x_date(date_labels = "%d-%b-%y", date_breaks = "1 month") +
  xlab("Date") +
  ylab("PM10 concentration") +
  labs(title = "Temporal Pattern of PM10 Concentration") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA)
  )


### Temporal Pattern of Air Pollutants ------------------------------------------
data_exclude_X <- pollution[, -1]

data_long <- melt(data_exclude_X, id.vars = "Date",
                  variable.name = "Pollutant",
                  value.name = "Concentration")

# Plot all variables in a single plot
ggplot(data_long, aes(x = Date, y = Concentration, color = Pollutant)) +
  geom_line() +
  scale_x_date(date_labels = "%d-%b-%y", date_breaks = "3 month") +
  labs(x = "Date", y = "Concentration", color = "Pollutant") +
  ggtitle("Temporal Pattern of Air Pollutants") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))


## Boxplots ---------------------------------------------------------------------
pm10.Box <- ggplot(pollution, aes(y = PM10)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("PM10 concentration")

carbonMonoxide.Box <- ggplot(pollution, aes(y = Carbon_Monoxide)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Carbon Monoxide concentration")

nitricOxigen.Box <- ggplot(pollution, aes(y = Nitric_Oxigen)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Nitric Oxide concentration")

nitrogenDioxide.Box <- ggplot(pollution, aes(y = Nitrogen_Dioxide)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Nitrogen Dioxide concentration")

BlackCarbon.Box <- ggplot(pollution, aes(y = BlackCarbon)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Black Carbon concentration")

Other.Box <- ggplot(pollution, aes(y = Other)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Other concentration")

Ozone.Box <- ggplot(pollution, aes(y = Ozone)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("Ozone concentration")

pm2.5.Box <- ggplot(pollution, aes(y = PM2.5)) +
  geom_boxplot(colour = "grey", fill = "palevioletred3") +
  ggtitle("PM2.5 concentration")


# Arrange the plots in a grid
grid.arrange(pm10.Box, carbonMonoxide.Box, nitricOxigen.Box, nitrogenDioxide.Box, BlackCarbon.Box,
             Other.Box, Ozone.Box, pm2.5.Box, ncol = 3)

## Correlation between variables -----------------------------------------------

cor_data <- pollution[, c("PM10", "BlackCarbon", "Carbon_Monoxide", "Nitric_Oxigen", 
                     "Nitrogen_Dioxide", "Other", "Ozone", "PM2.5")]

cor_matrix <- cor(cor_data)

ggcorrplot(cor_matrix,
           type = "lower",
           lab = TRUE,
           lab_size = 3.5,
           colors = c("red", "white", "#56B4E9"),
           title = "Correlogram of Dataset",
           ggtheme = theme_bw)

## Pairs Plot ------------------------------------------------------------------
custom_ggpairs <- function(pollution) {
  ggpairs(pollution,
          lower = list(continuous = "points", size = 0.5),
          upper = list(continuous = wrap("cor", digits = 2)),
          title = "Pairwise Correlations")
}

custom_ggpairs(as.data.frame(pollution[, 3:10]))


## VARIATION -------------------------------------------------------------------
pm10_sd <- sd(pollution$PM10)

pm10_var <- var(pollution$PM10)

summary(pollution$PM10)

hist <- ggplot(pollution, aes(x = PM10)) +
  geom_histogram(fill = "#41ab5d", color = "black") +
  labs(x = "PM10", y = "Frequency") +
  ggtitle("Distribution of PM10")

box <- ggplot(pollution, aes(y = PM10)) +
  geom_boxplot(fill = "#41ab5d", color = "black") +
  labs(y = "PM10") +
  ggtitle("Box Plot of PM10")

density <- ggplot(pollution, aes(x = PM10)) +
  geom_density(fill = "#41ab5d", color = "black") +
  labs(x = "PM10", y = "Density") +
  ggtitle("Density Plot of PM10")

grid.arrange(box, hist, density, ncol = 3)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# declare time-series variable

data <- ts(pollution$PM10, start = c(2014,1,5), frequency = 52)

# plotting time series object
autoplot(data)+ ggtitle("Temporal Pattern of PM10 Concentration") + 
  labs(x="Time", y="PM10 Concentration")

# ACF AND PACF ----------------------------------------------------

ggAcf(data) + ggtitle("ACF of PM10")
ggPacf(data) + ggtitle("Pacf of PM10")

# DIFFERENCE THE SERIES
dinf <- diff(data)

ggAcf(dinf) + ggtitle("ACF of PM10 (differenced) ")
ggPacf(dinf) + ggtitle("Pacf of PM10 (differenced)")

# Decompose the series

ts_decompose(data, type= "additive", showline=T)

# try STL decomposition
stl_decomp <- stl(data, s.window = "periodic")

plot(stl_decomp, col="#41ab5d")

## Non Stationary Tests ########################################################

# Perform the Dickey-Fuller test for trend -------------------------------------
adf.test(data) # indicating non- stationery

# Augmented Dickey-Fuller Test

# Dickey-Fuller = -2.9483, Lag order = 5, p-value = 0.1785
# alternative hypothesis: stationary

adf.test(dinf) # indicating stationery

# Augmented Dickey-Fuller Test

# Dickey-Fuller = -7.9811, Lag order = 5, p-value = 0.01
# alternative hypothesis: stationary

## Compute the Fourier transform ----------------------------------------------
fft_result <- fft(data)

# Plot the amplitude spectrum
plot(Mod(fft_result),
     type = "l",
     xlab = "Frequency",
     ylab = "Amplitude",
     main = "Amplitude Spectrum")

# Phillips Perron Test ---------------------------------------------------------
pp.test(data) #indicating stationary

# Phillips-Perron Unit Root Test

#Dickey-Fuller Z(alpha) = -79.584, Truncation lag parameter = 4, p-value = 0.01
# alternative hypothesis: stationary


# KPSS Test --------------------------------------------------------------------
kpss.test(data) # indicating stationary 

# KPSS Test for Level Stationarity
#KPSS Level = 0.087817, Truncation lag parameter = 4, p-value = 0.1


## In-sample Forecast ###########################################################

# Split the sample into Training and Test Sets

split_pm10 <- ts_split(data,sample= 52)

train_data <- split_pm10$train
test_data <- split_pm10$test

length(train_data)
length(test_data)

# Arima Diagnostic Plot on the Training set
arima_diag(train_data)

## removing trend and seasonality by differencing ------------------------------

train_data %>% diff() %>% ggtsdisplay(main="")

# - Plot new addmissions and detrended addmissions - #
par(mfrow = c(3,1))
plot(train_data, type = "l", col = "#41ab5d", main="Time Series") 
plot(diff(train_data, differences = 1), type = "l", col = "#41ab5d",main = "DIFF = 1") 
plot(diff(train_data, differences = 2), type = "l", col = "#41ab5d", main="DIFF = 2") 

## DIFF 1 --------------------------------------------------------------------

Diff.pm10 <- diff(train_data, differences = 1) # - difference the data
plot(Diff.pm10 , type = "l")

adf.test(Diff.pm10 , alternative = "stationary") # stationary

tsdisplay(Diff.pm10) # suggesting an (ma(3))


## remove trend and seasonality by linear + harmonic ---------------------------
x <- train_data
n <- length(x)
t <- 1:n

Z.fixed <- cbind(t, sin(2*pi*t/52), cos(2*pi*t/52))
trend.fixed <- lm(x~Z.fixed)$fitted.values
x.fixed <- x -trend.fixed

detrended <- ts(x.fixed, frequency = 52) 
tsdisplay(detrended, main="Detrended Time Series", col="#41ab5d") # seems like an AR(1) process but could also try ma(1).
adf.test(detrended)


### ARIMA Model Fits ###########################################################
# Auto ARIMA -------------------------------------------------------------------
fit.auto <- auto.arima(train_data, seasonal = T) # suggests order arima(0,1,2)(1,1,0)[52]
tsdisplay(residuals(fit.auto), col="#41ab5d", main="Arima(0,1,2)(1,1,0)")
summary(fit.auto)
check_res(fit.auto)
autoplot(fit.auto)


# Fit 1 -------------------------------------------------------------------------
fit1 <- arima(train_data, order = c(0,1,2), seasonal = c(1,1,0)) 
tsdisplay(residuals(fit1), lag.max = 55)
summary(fit1)
autoplot(fit1)
check_res(fit1)

# Fit 2 ------------------------------------------------------------------------
fit2 <- arima(train_data, order = c(0,1,2), seasonal = c(1,0,0)) 
tsdisplay(residuals(fit2), lag.max = 55)
summary(fit2)
autoplot(fit2)
check_res(fit2)

# Fit 3 ------------------------------------------------------------------------
fit3 <- arima(train_data, order = c(0,1,3), seasonal = c(1,1,0)) 
tsdisplay(residuals(fit3), lag.max = 55)
summary(fit3)
autoplot(fit3) 
check_res(fit3)

# Fit 4 Arima with harmonic regression -----------------------------------------
# # Set up harmonic regressors of order 2
fourier_terms <- 2

# Train the Auto ARIMA model with Fourier terms
fit4 <- auto.arima(train_data, seasonal = F , xreg = fourier(train_data, K = fourier_terms))
confint(fit4) # sin coeff are insignificant
forecast_periods <- length(test_data)
fcast4 <- forecast(fit4, h = forecast_periods, xreg = fourier(test_data, K = fourier_terms))

forecast_values <- as.vector(fcast4$mean)
lower_bound <- as.vector(fcast4$lower)
upper_bound <- as.vector(fcast4$upper)
accuracy_measures <- accuracy(fcast4, test_data)
print(accuracy_measures)


plot(forecast_result, main = "")
lines(test_data, col = "blue", lwd=1.5)
lines(fcast4$lower[, "80%"], col = "red", lty = 2)
lines(fcast4$upper[, "80%"], col = "red", lty = 2)
lines(fcast4$lower[, "95%"], col = "orange", lty = 2)
lines(fcast4$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

tsdisplay(residuals(fit4))
autoplot(fit4)

## Fit 5 with cosine only ------------------------------------------------------
# Create a sequence of 1 to the length of the training set (number of observations)
x <- train_data
t <- 1:length(train_data)

# Create two cosine functions with 52 weeks (frequency=1/52)
cosine_func_1 <- cos(2 * pi * t / 52)
cosine_func_2 <- cos(4 * pi * t / 52)

linear_model <- lm(train_data ~ t + cosine_func_1 + cosine_func_2)
summary(linear_model)
confint(linear_model)

# Train the Auto ARIMA model with the two cosine functions
xreg_data <- cbind(cosine_func_1, cosine_func_2)
fit5 <- auto.arima(train_data, seasonal= F, xreg = xreg_data)
summary(fit5)
confint(fit5)
coef(fit5)

tsdisplay(residuals(fit5), main="Arima Model")

forecast_periods <- length(test_data)
t_forecast <- (length(train_data) + 1):(length(train_data) + forecast_periods)
xreg_forecast <- cbind(cos(2 * pi * t_forecast / 52), cos(4 * pi * t_forecast / 52))
fcast5 <- forecast(fit5, h = forecast_periods, xreg = xreg_forecast)


forecast_values <- as.vector(fcast5$mean)
lower_bound <- as.vector(fcast5$lower)
upper_bound <- as.vector(fcast5$upper)
accuracy_measures <- accuracy(fcast5, test_data)

print(accuracy_measures)

# Extract the residuals from the ARIMA model 'fit5'
residuals_fit5 <- residuals(fit5)

ljung_box_test_fit5 <- Box.test(residuals_fit5, lag = 52, type = "Ljung-Box")

# Print the Ljung-Box test results
print(ljung_box_test_fit5)


plot(fcast5, main = "Arima Model")
lines(test_data, col = "blue", lwd=1.5)
lines(fcast5$lower[, "80%"], col = "red", lty = 2)
lines(fccast5$upper[, "80%"], col = "red", lty = 2)
lines(fcast5$lower[, "95%"], col = "orange", lty = 2)
lines(fcast5$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

autoplot(fit5)

## Comparing AIC ---------------------------------------------------------------
AIC(fit1, fit2, fit3, fit4, fit5)
BIC(fit1, fit2, fit3, fit4, fit5)

# Model        AIC
# fit1    847.8121
# fit2   1213.7417
# fit3    849.8088 
# fit4   1183.4850 
# fit5   1182.3451 

auto.fcast <- forecast(fit.auto, h=52)
fcast1 <- forecast(fit1, h=52)
fcast2 <- forecast(fit2, h=52)
fcast3 <- forecast(fit3, h=52)

## Accuracy --------------------------------------------------------------------

accuracy(auto.fcast, x = test_data)       # RMSE = 14.98
accuracy(fcast1, x = test_data)           # RMSE = 14.98
accuracy(fcast2, x = test_data)           # RMSE = 20.84
accuracy(fcast3, x = test_data)           # RMSE = 14.98
accuracy(fcast4, x = test_data)           # RMSE = 13.27
accuracy(fcast5, x = test_data)           # RMSE = 13.40

test_forecast(actual = data, forecast.obj = auto.fcast, test= test_data)
test_forecast(actual=data, forecast.obj = fcast1, test=test_data) 
test_forecast(actual=data, forecast.obj = fcast2, test=test_data)
test_forecast(actual=data, forecast.obj = fcast3, test=test_data)
test_forecast(actual=data, forecast.obj = fcast4, test=test_data)
test_forecast(actual=data, forecast.obj = fcast5, test= test_data)


# ACF and PACF plots for best fitting models ---------------------------------

p1<- autoplot(acf(fit.auto$residuals, plot=FALSE), main="ACF for ARIMA(0,1,2)(1,1,0)")
p2<- autoplot(pacf(fit.auto$residuals, plot=FALSE), main="PACF for ARIMA(0,1,2)(1,1,0)")
p3<- autoplot(pacf(fit5$residuals, plot=FALSE), main="ACF for Regression + ARIMA(0,1,2)")
p4<- autoplot(pacf(fit5$residuals, plot=FALSE), main="PACF for + Regression + ARIMA(0,1,2)")
grid.arrange(p1,p2,p3,p4, nrow=2)


### MODEL COMPARISON ###########################################################
## EXPONENTIAL SMOOTHING -------------------------------------------------------
## Holt-Winters Method 
holt_winters_model <- HoltWinters(train_data, seasonal = "additive")

forecast_periods <- length(test_data)

# Forecast using the trained model
forecast_result <- forecast(holt_winters_model, h = forecast_periods)

# Extract the forecasted values and prediction intervals
forecast_values <- as.vector(forecast_result$mean)
lower_bound <- as.vector(forecast_result$lower)
upper_bound <- as.vector(forecast_result$upper)

# Calculate accuracy measures
accuracy_measures <- accuracy(forecast_result, test_data)
accuracy_measures


# Plot
plot(train_data, xlim = c(2014, 2018),
     main = "Holt-Winters Fit", ylab = "PM10 Concentraion")
lines(hw.res.add$fitted[, 1], lty = 2, col = "blue")
lines(pred1.add[, 1], col = "blue", lwd = 1.5)
lines(pred1.add[, 2], col = "orange", lty = 2)
lines(pred1.add[, 3], col = "orange", lty = 2)
lines(test_data, col = "black", lty = 1)
legend(2014, 90,
       legend = c("Source", "HW Fit", "Prediction", "95% CI"),
       col = c("black", "blue", "blue", "orange"),
       lty = c(1, 2, 1, 2),
       lwd = c(1,1,3,1),
       cex = 0.5)

# Apply the checkresiduals() function
checkresiduals(holt_winters_model)


## Out of Sample Testing 36 weeks ##############################################
# Load the actual data for comparison
actual_data <- read.csv("http://www.stats.gla.ac.uk/~tereza/rp/Air_pollution_36.csv")
str(actual_data)
view(actual_data)

actual_data$Date <- as.Date(actual_data$Date,
                            format = "%Y-%m-%d %H:%M:%S")
data.add <- ts(actual_data$PM10, start = c(2018,1,7), frequency = 52)

## Auto Arima Out-of-Sample forecast -------------------------------------------
finalfit <- auto.arima(data)
autoplot(finalfit)
check_res(finalfit)
confint(finalfit)
coef(finalfit)

# Generate out of Sample Forecast
ffcast <- forecast(data, model=finalfit, h=36)
plot_forecast(ffcast)
accuracy(ffcast, x = data.add)

plot(ffcast, main = "Auto Arima Model- 36 weeks ahead Forecast")
lines(data.add, col = "blue", lwd=1.5)
lines(ffcast$lower[, "80%"], col = "red", lty = 2)
lines(ffcast$upper[, "80%"], col = "red", lty = 2)
lines(ffcast$lower[, "95%"], col = "orange", lty = 2)
lines(ffcast$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

## Final model Out-of-Sample Forecast ------------------------------------------
t <- 1:length(data)

# Create two cosine functions with 52 weeks (frequency=1/52)
cosine_func_1 <- cos(2 * pi * t / 52)
cosine_func_2 <- cos(4 * pi * t / 52)

# Train the Auto ARIMA model with the two cosine functions
xreg_data <- cbind(cosine_func_1, cosine_func_2)
arima_model <- auto.arima(data, seasonal = F, xreg = xreg_data)

forecast_periods <- 36
t_forecast <- (length(data) + 1):(length(data) + forecast_periods)
xreg_forecast <- cbind(cos(2 * pi * t_forecast / 52), cos(4 * pi * t_forecast / 52))
forecast_result <- forecast(arima_model, h = forecast_periods, xreg = xreg_forecast)


forecast_values <- as.vector(forecast_result$mean)
lower_bound <- as.vector(forecast_result$lower)
upper_bound <- as.vector(forecast_result$upper)
accuracy_measures <- accuracy(forecast_result, data.add)
print(accuracy_measures)

plot(forecast_result, main = "Regression with Arima Model - 36 weeks ahead Forecast")
lines(data.add, col = "blue", lwd = 1.5)
lines(forecast_result$lower[, "80%"], col = "red", lty = 2)
lines(forecast_result$upper[, "80%"], col = "red", lty = 2)
lines(forecast_result$lower[, "95%"], col = "orange", lty = 2)
lines(forecast_result$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

plot_forecast(forecast_result)

## Out-of-Sample Forecast 120 weeks ############################################
t <- 1:length(data)

# Create two cosine functions with 52 weeks (frequency=1/52)
cosine_func_1 <- cos(2 * pi * t / 52)
cosine_func_2 <- cos(4 * pi * t / 52)

# Train the Auto ARIMA model with the two cosine functions
xreg_data <- cbind(cosine_func_1, cosine_func_2)
arima_model <- auto.arima(data, seasonal = F, xreg = xreg_data)

forecast_periods <- 120
t_forecast <- (length(data) + 1):(length(data) + forecast_periods)
xreg_forecast <- cbind(cos(2 * pi * t_forecast / 52), cos(4 * pi * t_forecast / 52))
forecast_result <- forecast(arima_model, h = forecast_periods, xreg = xreg_forecast)


forecast_values <- as.vector(forecast_result$mean)
lower_bound <- as.vector(forecast_result$lower)
upper_bound <- as.vector(forecast_result$upper)
accuracy_measures <- accuracy(forecast_result, data.add)
print(accuracy_measures)

plot(forecast_result, main = "Arima Model - 120 weeks ahead Forecast")
lines(data.add, col = "blue", lwd = 1.5)
lines(forecast_result$lower[, "80%"], col = "red", lty = 2)
lines(forecast_result$upper[, "80%"], col = "red", lty = 2)
lines(forecast_result$lower[, "95%"], col = "orange", lty = 2)
lines(forecast_result$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

### Auto arima 120 weeks ahead forecast ----------------------------------------
ffcast <- forecast(data, model=finalfit, h=120)
plot_forecast(ffcast)
accuracy(ffcast, x = data.add)

plot(ffcast, main = "Auto Arima Model- 120 weeks ahead Forecast")
lines(data.add, col = "blue", lwd=1.5)
lines(ffcast$lower[, "80%"], col = "red", lty = 2)
lines(ffcast$upper[, "80%"], col = "red", lty = 2)
lines(ffcast$lower[, "95%"], col = "orange", lty = 2)
lines(ffcast$upper[, "95%"], col = "orange", lty = 2)
legend("topleft", legend = c("Series", "80% CI", "95% CI", "Actual"),
       col = c("black", "orange", "red", "blue"), 
       lty = c(1, 2, 2, 1),
       lwd = c(1, 1, 1, 2),
       cex = 0.5)

