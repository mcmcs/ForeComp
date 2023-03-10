

Prepare_SPF_Data <- function(series_name, start_year, end_year) {

  df_gdp <- read_csv(str_glue("{series_name}_extended.csv")) %>%
    select(X1, starts_with("SPF"), starts_with("NC"), Realiz1) %>%
    mutate(
      year_quarter  = str_replace(X1, ":", "-") %>% yq(.),
      across(-year_quarter, ~ as.double(.))
    ) %>%
    select(-X1) %>%
    select(year_quarter, everything()) %>%
    filter(between(year_quarter, as_date(start_year), as_date(end_year))) %>%
    glimpse()

  e1 <- df_gdp$Realiz1 - df_gdp$NCfor_Step1
  e1[is.na(e1)] <- 0

  e2 <- df_gdp$Realiz1 - df_gdp$SPFfor_Step1
  e2[is.na(e2)] <- 0

  d <- (e1 ^ 2) - (e2 ^ 2)

}
