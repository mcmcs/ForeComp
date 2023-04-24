
pacman::p_load(
  tidyverse,
  lubridate,
  tidylog,
  forecast,
  astsa,
  furrr,
  tictoc
)

library(ForeComp)


source("Prepare_SPF_Data.R")
source("Plot_Power_Size_Tradeoff.R")

v_series <- c("RGDP", "TBILL", "UNEMP", "PGDP")
v_horizon <- 1:5

df_start_end <- tibble(
  starting = c("1987-01-01", "1987-01-01", "1987-01-01", "1997-01-01", "2007-01-01", "2017-01-01"),
  ending = c("2016-12-01", "2021-12-01", "1996-12-01", "2006-12-01", "2016-12-01", "2021-12-01")
)

df_combo <- expand_grid(
  series = v_series,
  horizon = v_horizon,
  df_start_end
)

set.seed(1234)

for (combo_index in 1:nrow(df_combo)) {

  series <- df_combo %>%
    slice(combo_index) %>%
    pull(series)

  starting_year <- df_combo %>%
    slice(combo_index) %>%
    pull(starting)

  ending_year <- df_combo %>%
    slice(combo_index) %>%
    pull(ending)

  horizon <- df_combo %>%
    slice(combo_index) %>%
    pull(horizon)

  print(str_glue("Running calcations for row {combo_index}/{nrow(df_combo)}"))

  df_data_prepped <- Prepare_SPF_Data(series, horizon, starting_year, ending_year)
  Plot_Size_Power_Tradeoff(raw_data = df_data_prepped,
                           nlen = length(df_data_prepped), nsim = 10000,
                           cl = .05,
                           M_set = c(1:10, seq(11, length(df_data_prepped) - 1, 10)))
}






