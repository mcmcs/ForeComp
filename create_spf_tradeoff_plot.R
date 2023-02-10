
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

nsim_ <- 10000
cl_ <- .05

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

mu_t <- mean(d)
v_dt_tilde <- d - mu_t

stopifnot(all.equal(0, mean(v_dt_tilde)))

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

del_grid <- seq(0, 4, .05)

# Size Computation --------------------------------------------------------

l_arima_sim <- map(1:nsim_, ~ arima.sim(m_sim, n = n, innov = rnorm(n, 0, sqrt(a$sigma2))))

v_oracle_test_statistic <- map_dbl(l_arima_sim, ~ mean(.) / sqrt(Om / n))
v_oracle_p_val <-  2 * stats::pnorm(-abs(v_oracle_test_statistic), mean=0, sd=1)
v_oracle_rej <- v_oracle_p_val < cl_
empirical_size_oracle_ <- mean(v_oracle_rej)
c05_star_oracle <- quantile(abs(v_oracle_test_statistic), 1 - cl_)

# size corrected power ----------------------------------------------------

df_arima_del_grid <- expand_grid(l_arima_sim, del_grid) %>%
  mutate(
    oracle_size_corrected_stat = map2(l_arima_sim, del_grid, ~ .x + (1/sqrt(n)) *sqrt(Om)* .y),
    oracle_dm_stat = map_dbl(oracle_size_corrected_stat, ~ mean(.) / (sqrt(Om / n))),
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

# 87 seconds for 14 bandwidths w/ 1000 sims

v_M <- c(1:10, seq(20, 50, 10))

plan(multisession, workers = 4)

options(future.globals.maxSize= 10 * 1073741824)

tic()
df_sims <- future_map_dfr(v_M, ~ Compute_Size_Distort_Max_Power_Loss(l_arima_sim, .)) %>%
  print()
toc()

ggplot(df_sims, aes(x = size_distortion, y = max_power_loss, label = M)) +
  geom_point(size = 6) +
  geom_text(nudge_x = .003, nudge_y = .001) +
  theme_minimal() +
  labs(
    x = "Size Distortion",
    y = "Max Power Loss"
  ) +
  theme(
    panel.grid = element_blank()
  )



