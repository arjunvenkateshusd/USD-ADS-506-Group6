# Libraries we use
library(tidyverse)   # data manipulation / plotting
library(lubridate)   # date time stuff
library(forecast)    # time series tools
library(zoo)         # rolling means
library(tseries)     # adf test


# load up our dataset
data_raw <- read.csv("processed_ems_load_data.csv", stringsAsFactors = FALSE)

str(data_raw)
summary(data_raw)

# parse timestamps
data_raw <- data_raw %>%
  mutate(timestamp = ymd_hms(timestamp))

data <- data_raw %>%
  select(
    timestamp,
    HR,
    PGE, SCE, SDGE, VEA, CAISO,
    hour, day, weekday, month, year, weekofyear, is_weekend, season
  )

range(data$timestamp)
nrow(data)
n_distinct(data$timestamp)

# eda
# 1. missing values per column
data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  print()

# 2. stats for load by region
data %>%
  summarise(
    PGE_min = min(PGE),  PGE_mean = mean(PGE),  PGE_max = max(PGE),
    SCE_min = min(SCE),  SCE_mean = mean(SCE),  SCE_max = max(SCE),
    SDGE_min = min(SDGE), SDGE_mean = mean(SDGE), SDGE_max = max(SDGE),
    VEA_min = min(VEA),  VEA_mean = mean(VEA),  VEA_max = max(VEA),
    CAISO_min = min(CAISO), CAISO_mean = mean(CAISO), CAISO_max = max(CAISO)
  ) %>%
  print()

# 3. for multi region plots
load_long <- data %>%
  select(timestamp, PGE, SCE, SDGE, VEA, CAISO) %>%
  pivot_longer(
    cols = -timestamp,
    names_to = "region",
    values_to = "load"
  )

# 4. hourly load by region
ggplot(load_long, aes(x = timestamp, y = load)) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  labs(
    title = "Hourly Electricity Load by Region",
    x = "Time",
    y = "Load (MW)"
  )

# 5. daily average load by region
daily_load <- load_long %>%
  mutate(date = as.Date(timestamp)) %>%
  group_by(region, date) %>%
  summarise(mean_load = mean(load, na.rm = TRUE), .groups = "drop")

ggplot(daily_load, aes(x = date, y = mean_load)) +
  geom_line() +
  facet_wrap(~ region, scales = "free_y", ncol = 2) +
  labs(
    title = "Daily Average Load by Region",
    x = "Date",
    y = "Average Load (MW)"
  )

# 6. hour of day load pattern (CAISO)
ggplot(
  data %>% mutate(hour_factor = factor(hour)),
  aes(x = hour_factor, y = CAISO)
) +
  geom_boxplot() +
  labs(
    title = "CAISO Load by Hour of Day",
    x = "Hour of Day",
    y = "Load (MW)"
  )

# 7. seasonal pattern for CAISO
ggplot(
  data %>% mutate(season = factor(season)),
  aes(x = hour, y = CAISO, group = interaction(season, day), color = season)
) +
  geom_line(alpha = 0.1) +
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = season),
    linewidth = 1.1
  ) +
  labs(
    title = "CAISO Hourly Load by Season (Mean Overlay)",
    x = "Hour of Day",
    y = "Load (MW)"
  )


# preprocessing for time series

data_ts <- data %>%
  arrange(timestamp)

# hourly CAISO series (24 obs per day)
ts_caiso <- ts(data_ts$CAISO, frequency = 24)

# 1. raw CAISO time series
autoplot(ts_caiso) +
  labs(
    title = "CAISO Hourly Load (Raw Series)",
    x = "Time (hours)",
    y = "Load (MW)"
  )

# 2. stl decomposition
caiso_stl <- stl(ts_caiso, s.window = "periodic")
autoplot(caiso_stl)

# 3. smoothing via 24 hour moving average
caiso_ma24 <- ma(ts_caiso, order = 24)

autoplot(caiso_ma24) +
  labs(
    title = "CAISO Smoothed Load (24-hour Moving Average)",
    x = "Time (hours)",
    y = "Smoothed Load (MW)"
  )

# 4. box-cox
lambda_caiso <- BoxCox.lambda(ts_caiso)
lambda_caiso

ts_caiso_bc <- BoxCox(ts_caiso, lambda_caiso)

autoplot(ts_caiso_bc) +
  labs(
    title = "CAISO Box-Cox Transformed Load",
    x = "Time (hours)",
    y = "Transformed Load"
  )

# 5. differencing for stationarity

ts_caiso_bc_diff1 <- diff(ts_caiso_bc, lag = 1)

autoplot(ts_caiso_bc_diff1) +
  labs(
    title = "CAISO Box-Cox Transformed Load (First Difference)",
    x = "Time (hours)",
    y = "Differenced Load"
  )

# seasonal difference
ts_caiso_bc_diff1_seas <- diff(ts_caiso_bc_diff1, lag = 24)

autoplot(ts_caiso_bc_diff1_seas) +
  labs(
    title = "CAISO Box-Cox Transformed Load (First + Seasonal Difference)",
    x = "Time (hours)",
    y = "Differenced Load"
  )

# 6. acf/pacf
ts_caiso_clean <- ts_caiso_bc_diff1_seas[!is.na(ts_caiso_bc_diff1_seas)]

ggAcf(ts_caiso_clean) +
  labs(title = "ACF: CAISO Differenced Series")

ggPacf(ts_caiso_clean) +
  labs(title = "PACF: CAISO Differenced Series")

adf_result <- adf.test(ts_caiso_clean)
adf_result


# df

model_df <- data_ts %>%
  mutate(
    CAISO_bc = as.numeric(ts_caiso_bc),
    CAISO_bc_diff1 = c(NA, as.numeric(ts_caiso_bc_diff1)),
    CAISO_bc_diff1_seas = c(rep(NA, 1 + 24), as.numeric(ts_caiso_clean))
  )

head(model_df)

# if you want to uncomment the line below for caiso_preprocessed_for_modeling.csv
# write.csv(model_df, "caiso_preprocessed_for_modeling.csv", row.names = FALSE)



# Arjun's code above - below I will create the ARIMA 

data_raw <- read.csv("processed_ems_load_data.csv")

data_raw <- data_raw %>%
  mutate(timestamp = ymd_hms(timestamp)) %>%
  arrange(timestamp)

# to 24 hrs
caiso_ts <- ts(data_raw$CAISO, frequency = 24)

#train / test 
n <- length(caiso_ts)
n_train <- floor(0.8 * n)

train_ts <- ts(caiso_ts[1:n_train], frequency = 24)
test_ts  <- ts(caiso_ts[(n_train + 1):n], frequency = 24)

length(train_ts)
length(test_ts)

#auto arima 
fit_auto <- auto.arima(train_ts,
                       seasonal = TRUE,
                       approximation = TRUE)
fit_auto
fc_auto <- forecast(fit_auto, h = length(test_ts))

accuracy(fc_auto, test_ts)

#manual arima 
fit_manual <- Arima(train_ts,
                    order = c(1,1,1),
                    seasonal = c(1,1,1))

fit_manual
fc_manual <- forecast(fit_manual, h = length(test_ts))

accuracy(fc_manual, test_ts)

#plot 
autoplot(fc_auto) +
  autolayer(test_ts, series = "Test")

autoplot(fc_manual) +
  autolayer(test_ts, series = "Test")
