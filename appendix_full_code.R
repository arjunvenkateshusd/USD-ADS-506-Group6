# ============================
# Libraries we use
# ============================
library(tidyverse)   # data manipulation / plotting
library(lubridate)   # date time stuff
library(forecast)    # time series tools
library(zoo)         # rolling means
library(tseries)     # ADF test
library(kableExtra)  # APA-style tables
library(patchwork)   # combine ggplots

# ============================
# Load dataset
# ============================
data_raw <- read.csv("processed_ems_load_data.csv", stringsAsFactors = FALSE)

str(data_raw)
summary(data_raw)

# Parse timestamps
data_raw <- data_raw %>%
  mutate(timestamp = ymd_hms(timestamp))

data <- data_raw %>%
  select(
    timestamp, HR, PGE, SCE, SDGE, VEA, CAISO,
    hour, day, weekday, month, year, weekofyear, is_weekend, season
  )

range(data$timestamp)
nrow(data)
n_distinct(data$timestamp)

# ============================
# Exploratory Data Analysis
# ============================

# 1. Missing values per column
data %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  print()

# 2. Summary statistics for regions
data %>%
  summarise(
    PGE_min = min(PGE),  PGE_mean = mean(PGE),  PGE_max = max(PGE),
    SCE_min = min(SCE),  SCE_mean = mean(SCE),  SCE_max = max(SCE),
    SDGE_min = min(SDGE), SDGE_mean = mean(SDGE), SDGE_max = max(SDGE),
    VEA_min = min(VEA),  VEA_mean = mean(VEA),  VEA_max = max(VEA),
    CAISO_min = min(CAISO), CAISO_mean = mean(CAISO), CAISO_max = max(CAISO)
  ) %>% print()

# 3. Reshape for multi-region EDA plots
load_long <- data %>%
  select(timestamp, PGE, SCE, SDGE, VEA, CAISO) %>%
  pivot_longer(cols = -timestamp, names_to = "region", values_to = "load")

# 4. Daily average load by region
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

# 5. CAISO hourly load pattern
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

# 6. Seasonal CAISO pattern
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

# ============================
# Time Series Object + Preprocessing
# ============================

data_ts <- data %>% arrange(timestamp)
ts_caiso <- ts(data_ts$CAISO, frequency = 24)

# STL decomposition (trend + seasonality)
caiso_stl <- stl(ts_caiso, s.window = "periodic")
autoplot(caiso_stl)

# 24-hour moving average smoothing
caiso_ma24 <- ma(ts_caiso, order = 24)
autoplot(caiso_ma24) +
  labs(
    title = "CAISO Smoothed Load (24-Hour Moving Average)",
    x = "Time (hours)",
    y = "Smoothed Load (MW)"
  )

# Box-Cox transformation
lambda_caiso <- BoxCox.lambda(ts_caiso)
ts_caiso_bc <- BoxCox(ts_caiso, lambda_caiso)

autoplot(ts_caiso_bc) +
  labs(
    title = "CAISO Box-Cox Transformed Load",
    x = "Time (hours)",
    y = "Transformed Load"
  )

# Differencing for stationarity: first + seasonal
ts_caiso_bc_diff1 <- diff(ts_caiso_bc, lag = 1)

autoplot(ts_caiso_bc_diff1) +
  labs(
    title = "CAISO Box-Cox Transformed Load (First Difference)",
    x = "Time (hours)",
    y = "Differenced Load"
  )

ts_caiso_bc_diff1_seas <- diff(ts_caiso_bc_diff1, lag = 24)

autoplot(ts_caiso_bc_diff1_seas) +
  labs(
    title = "CAISO Box-Cox Transformed Load (First + Seasonal Difference)",
    x = "Time (hours)",
    y = "Differenced Load"
  )

ts_caiso_clean <- ts_caiso_bc_diff1_seas[!is.na(ts_caiso_bc_diff1_seas)]

# ACF/PACF on differenced series
ggAcf(ts_caiso_clean) + labs(title = "ACF: Differenced CAISO Series")
ggPacf(ts_caiso_clean) + labs(title = "PACF: Differenced CAISO Series")

# ADF test for stationarity
adf_result <- adf.test(ts_caiso_clean)
adf_result

# ============================
# Train–Test Split
# ============================
n <- length(ts_caiso)
n_train <- floor(0.8 * n)

train_ts <- ts(ts_caiso[1:n_train], frequency = 24)
test_ts  <- ts(ts_caiso[(n_train + 1):n], frequency = 24)

length(train_ts)
length(test_ts)

# ============================
# Auto ARIMA
# ============================
fit_auto <- auto.arima(train_ts,
                       seasonal = TRUE,
                       approximation = TRUE)

fit_auto
fc_auto <- forecast(fit_auto, h = length(test_ts))

# Use numeric test vector to avoid ts window() issues
acc_auto <- accuracy(fc_auto, as.numeric(test_ts))
acc_auto

# ============================
# Manual ARIMA (fixed to be stable and comparable)
# ============================
# Using structure inspired by auto.arima and ACF/PACF:
# ARIMA(5,0,0)(2,1,0)[24]

fit_manual <- Arima(
  train_ts,
  order = c(2, 1, 2),
  seasonal = list(order = c(1, 1, 1), period = 24)
)

fit_manual
fc_manual <- forecast(fit_manual, h = length(test_ts))

acc_manual <- accuracy(fc_manual, as.numeric(test_ts))
acc_manual

# ============================
# TSLM model
# ============================
fit_tslm <- tslm(train_ts ~ trend + season)
summary(fit_tslm)

fc_tslm <- forecast(fit_tslm, h = length(test_ts))
acc_tslm <- accuracy(fc_tslm, as.numeric(test_ts))
acc_tslm

# ============================
# ETS model
# ============================
fit_ets <- ets(train_ts)
fit_ets

fc_ets <- forecast(fit_ets, h = length(test_ts))
acc_ets <- accuracy(fc_ets, as.numeric(test_ts))
acc_ets

# ============================
# Forecast comparison plot
# ============================
autoplot(fc_auto) +
  autolayer(fc_manual$mean, series = "ARIMA (Manual)") +
  autolayer(fc_tslm$mean,   series = "TSLM (Trend + Season)") +
  autolayer(fc_ets$mean,    series = "ETS") +
  autolayer(test_ts,        series = "Actual Load") +
  labs(
    title = "Forecast Comparison: ARIMA, TSLM, ETS vs. Actual CAISO Load",
    x = "Time",
    y = "Load (MW)"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 16, face = "bold"),
    legend.position = "bottom"
  )

# ============================
# APA 7 Model Accuracy Table
# ============================
accuracy_table <- rbind(
  `ARIMA (Auto)`   = acc_auto["Test set", ],
  `ARIMA (Manual)` = acc_manual["Test set", ],
  `ETS`            = acc_ets["Test set", ],
  `TSLM`           = acc_tslm["Test set", ]
)

accuracy_table %>%
  kable(
    format  = "html",
    caption = "Table 1\n\nModel Accuracy Comparison for CAISO Load Forecasting",
    digits  = 3
  ) %>%
  kable_styling(full_width = FALSE,
                position   = "center",
                font_size  = 12) %>%
  row_spec(0, bold = TRUE)

# ============================
# Side-by-side ACF + PACF using patchwork
# ============================
p1 <- ggAcf(ts_caiso_clean) +
  labs(title = "ACF After Seasonal Differencing")

p2 <- ggPacf(ts_caiso_clean) +
  labs(title = "PACF After Seasonal Differencing")

p1 + p2

    