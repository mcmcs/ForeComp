
pacman::p_load(
  tidyverse,
  lubridate
)

UNEMP <- read_csv("spf-data-raw/UNEMP_extended.csv") %>%
  select(1, starts_with("SPF"), starts_with("NC"), starts_with("Realiz")) %>%
  mutate(
    year_quarter = str_replace(...1, ":", "-") %>% yq(.),
    across(-year_quarter, ~ as.double(.)),
    across(-year_quarter, ~ coalesce(., 0))
  ) %>%
  select(-1) %>%
  select(year_quarter, everything())

usethis::use_data(UNEMP, overwrite = TRUE)
