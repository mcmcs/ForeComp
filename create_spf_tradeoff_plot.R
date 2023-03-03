
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

df_gdp <- read_csv("PGDP_extended.csv") %>%
    select(X1, starts_with("SPF"), starts_with("NC"), Realiz1) %>%
    mutate(
      year_quarter  = str_replace(X1, ":", "-") %>% yq(.),
      across(-year_quarter, ~ as.double(.))
    ) %>%
    select(-X1) %>%
    select(year_quarter, everything()) %>%
    filter(between(year_quarter, as_date("1987-01-01"), as_date("2016-10-01"))) %>%
    glimpse()

e1 <- df_gdp$Realiz1 - df_gdp$NCfor_Step1
e1[is.na(e1)] <- 0

e2 <- df_gdp$Realiz1 - df_gdp$SPFfor_Step1
e2[is.na(e2)] <- 0

d <- (e1 ^ 2) - (e2 ^ 2)
n <- length(d)

source("Plot_Power_Size_Tradeoff.R")

Plot_Size_Power_Tradeoff(raw_data = d,
                         nlen = n, nsim = 10000,
                         cl = .05,
                         M_set = c(seq(from=2, to=10, by=1), seq(from=11, to=150, by=10))[1:19])


