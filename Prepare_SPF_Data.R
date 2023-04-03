

Prepare_SPF_Data <- function(series, horizon, start_year, end_year) {

  spf_var <- str_glue("SPFfor_Step{horizon}")
  nc_var <- str_glue("NCfor_Step{horizon}")
  realiz_var <- str_glue("Realiz{horizon}")

  series_name <- str_glue("{series}_extended.csv")

  df_gdp <- read_csv(series_name) %>%
    select(X1, all_of(c(spf_var, nc_var, realiz_var))) %>%
    mutate(
      year_quarter  = str_replace(X1, ":", "-") %>% yq(.),
      across(-year_quarter, ~ as.double(.))
    ) %>%
    select(-X1) %>%
    filter(between(year_quarter, as_date(start_year), as_date(end_year))) %>%
    glimpse()

  e1 <- df_gdp[[realiz_var]] - df_gdp[[nc_var]]
  e1[is.na(e1)] <- 0

  e2 <- df_gdp[[realiz_var]] - df_gdp[[spf_var]]
  e2[is.na(e2)] <- 0

  d <- (e1 ^ 2) - (e2 ^ 2)

}
