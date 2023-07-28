
pacman::p_load(
  tidyverse,
  lubridate
)

TBILL <- read_csv("spf-data-raw/TBILL_extended.csv") %>%
  select(1, starts_with("SPF"), starts_with("NC"), starts_with("Realiz")) %>%
  mutate(
    year_quarter = str_replace(...1, ":", "-") %>% yq(.),
    across(-year_quarter, ~ as.double(.))
  ) %>%
  select(-1) %>%
  select(year_quarter, everything())

usethis::use_data(TBILL, overwrite = TRUE)
