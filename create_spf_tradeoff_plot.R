
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

df_combo <- tibble(
  series = v_series,
  starting = "1987-01-01",
  ending = "2020-12-01"
)

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

  print(str_glue("Running calcations for {series}"))

  df_data_prepped <- Prepare_SPF_Data(series, starting_year, ending_year)
  Plot_Size_Power_Tradeoff(raw_data = df_data_prepped,
                           nlen = length(df_data_prepped), nsim = 10000,
                           cl = .05,
                           M_set = c(seq(from=2, to=10, by=1), seq(from=11, to=150, by=10))[1:19])
}






