
pacman::p_load(
  tidyverse,
  forecast,
  astsa,
  ForeComp,
  tictoc,
  furrr
)

# Global Tuning Parameters ------------------------------------------------

nlen_ <- 200
nsim_ <- 10000
cl_ <- .05
ar1 <- .7

b <- seq(.1, 1, .1)
m <- floor(b * nlen_)

v_Mchoice <- c(1:9, m)

set.seed(1234)


# Simulated Data ----------------------------------------------------------

T_ <- 1000 + nlen_

dt0_ <- 0
v_dt <- vector("double", T_)

for (t in 1:T_) {

  dt0_ <- .5 + ar1 * dt0_ + rnorm(1)
  v_dt[t] <- dt0_

}

v_dt <- v_dt[1001:T_]

mu_t_ <- mean(v_dt)
v_dt_tilde <- v_dt - mu_t_

# ARIMA -------------------------------------------------------------------

a = auto.arima(y=v_dt_tilde, max.p = 12, max.q=12, stationary=T, ic="aic", seasonal=F, allowmean=F);
i_ar = grep("ar", names(a$coef));
i_ma = grep("ma", names(a$coef));
m_sim = list("ar"=a$coef[i_ar], "ma"=a$coef[i_ma]);
if (length(m_sim$ar)==0){m_sim$ar=0.0}
if (length(m_sim$ma)==0){m_sim$ma=0.0}

ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100);
Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than a spectral density at zero

# Delta setup -------------------------------------------------------------

ndel <-  50 #number of deltas
del_tilde <-  10 #largest delta
del_grid <-  seq(from=-del_tilde, to=del_tilde, length.out=ndel)

# Size Computation --------------------------------------------------------

l_arima_sim <- map(1:nsim_, ~ arima.sim(m_sim, n = nlen_, innov = rnorm(nlen_, 0, sqrt(a$sigma2))))

v_oracle_test_statistic <- map_dbl(l_arima_sim, ~ mean(.) / sqrt(Om / nlen_))
v_oracle_p_val <-  2 * stats::pnorm(-abs(v_oracle_test_statistic), mean=0, sd=1)
v_oracle_rej <- v_oracle_p_val < cl_
empirical_size_oracle_ <- mean(v_oracle_rej)
c05_star_oracle <- quantile(abs(v_oracle_test_statistic), 1 - cl_)

# size corrected power ----------------------------------------------------

df_arima_del_grid <- expand_grid(l_arima_sim, del_grid) %>%
  mutate(
    oracle_size_corrected_stat = map2(l_arima_sim, del_grid, ~ .x + (1/sqrt(nlen_)) *sqrt(Om)* .y),
    oracle_dm_stat = map_dbl(oracle_size_corrected_stat, ~ mean(.) / (sqrt(Om / nlen_))),
    oracle_pval = 2*stats::pnorm(-abs(oracle_dm_stat), mean=0, sd=1),
    oracle_reject_raw_power = oracle_pval < cl_,
    oracle_reject_size_corrected_power = abs(oracle_dm_stat) > c05_star_oracle
  )

Compute_Size_Distort_Max_Power_Loss <- function(simulated_data, bandwidth) {

  print(paste("Running simulation for M =", bandwidth))

  l_dm_test <- map(simulated_data, ~ dm.test.bt(., M = bandwidth, cl = cl_))
  v_test_statistic <- map_dbl(l_dm_test, ~ pluck(., "stat"))
  v_reject <- map_dbl(l_dm_test, ~ pluck(., "rej"))

  size_distortion_ <- mean(v_reject) - cl_

  v_standard_stat <- map(l_arima_sim, ~ dm.test.bt(., M = bandwidth, cl = cl_)) %>%
    map_dbl(., ~ pluck(., "stat"))

  c05_star <- quantile(abs(v_standard_stat), (1 - cl_))

  df_arima_del_grid <- df_arima_del_grid %>%
    mutate(
      standard_dm_test =  map(oracle_size_corrected_stat, ~ dm.test.bt(., M = bandwidth, cl = cl_)),
      standard_raw_power = map_lgl(standard_dm_test, ~ pluck(., "rej")),
      standard_size_corrected_stat = map_dbl(standard_dm_test, ~ pluck(., "stat")),
      standard_size_corrected_power = abs(standard_size_corrected_stat) > c05_star
    )

  df_power_loss <- df_arima_del_grid %>%
    group_by(del_grid) %>%
    summarize(
      standard_power = mean(standard_size_corrected_power),
      oracle_power = mean(oracle_reject_size_corrected_power),
      power_loss = oracle_power - standard_power
    )

  max_power_loss_ <- max(df_power_loss$oracle_power - df_power_loss$standard_power)

  tibble(
    M = bandwidth,
    size_distortion = size_distortion_,
    max_power_loss = max_power_loss_
  )



}

# 30 minutes

v_M <- 1:10

plan(multisession, workers = 5)

options(future.globals.maxSize= 10 * 1073741824)

tic()
df_sims <- future_map_dfr(v_M, ~ Compute_Size_Distort_Max_Power_Loss(l_arima_sim, .)) %>%
  print()
toc()

hyper_parameters <- pas

ggplot(filter(df_sims), aes(x = size_distortion, y = max_power_loss, label = M)) +
  geom_point(size = 6) +
  # geom_text(nudge_x = .003) +
  theme_minimal() +
  labs(
    x = "Size Distortion",
    y = "Max Power Loss"
  ) +
  theme(
    panel.grid = element_blank()
  )

file_out_ <- paste("data/sim", ar1, nlen_, nsim_, sep = "_")
write_csv(df_sims, paste0(file_out_, ".csv"))


# Max Power Loss ----------------------------------------------------------

 %>%
  ungroup() %>%
  mutate(
    d_max_power_loss = power_loss == max(power_loss)
  ) %>%
  print()

max_power_loss_ <-



df_power_loss %>%
  select(del_grid, standard_power, oracle_power) %>%
  pivot_longer(., -del_grid, names_to = "test_type", values_to = "power") %>%
  mutate(
    `Test Type` = str_remove(test_type, "_power"),
    `Test Type` = str_to_title(`Test Type`)
  ) %>%
  ggplot(., aes(x = del_grid, y = power)) +
  ggthemes::scale_color_ptol() +
  geom_line(aes(color = `Test Type`), size = 1) +
  geom_point(aes(color = `Test Type`), size = 1.5) +
  geom_segment(data = df_power_loss %>% filter(d_max_power_loss),
               aes(x = del_grid, xend = del_grid, y = standard_power, yend = oracle_power), size = 2) +
  theme_minimal() +
  labs(
    x = "Delta Grid",
    y = "Power"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12)
  )

paste("size distortion =", size_distortion_)


