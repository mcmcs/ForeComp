
pacman::p_load(tidyverse)

NUM_PERIODS <- c(8, 128)
HORIZON <- c(1, 3)
DM_TEST <- c("original", "modified", "ewc", "bk")

set.seed(1234)

Take_Draw <- function(num_periods, horizon, dm_test) {

  e1 <- rnorm(n = num_periods)
  e2 <- rnorm(n = num_periods)

  d = e1^2 - e2^2

  hypothesis_result <- switch(dm_test,
                    "original" = dm.test.r(d, h = horizon, cl = .10)$rej,
                    "modified" = dm.test.r.m(d, h = horizon, cl = .10)$rej,
                    "ewc" = dm.test.ewc.fb(d, cl = .10)$rej,
                    "bk" = dm.test.bt.fb(d, cl = .10)$rej
  )

  return(hypothesis_result)

}

(df_experiments <- expand_grid(
  num_periods = NUM_PERIODS,
  horizon = HORIZON,
  dm_test = DM_TEST)
)

l_simulations <- pmap(df_experiments, ~ rerun(10000, Take_Draw(..1, ..2, ..3)))

Calc_Reject_Rate <- function(simulated_values, index) {

  simulated_values %>%
    pluck(index) %>%
    unlist() %>%
    mean(., na.rm = TRUE) %>%
    `*`(., 100) %>%
    round(1)
}

(df_experiments <- df_experiments %>%
  mutate(
    reject_rate = map_dbl(1:nrow(df_experiments), ~ Calc_Reject_Rate(l_simulations, .))
  )
)

df_experiments %>%
  pivot_wider(., names_from = dm_test, values_from = reject_rate)

