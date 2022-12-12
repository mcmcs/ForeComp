
pacman::p_load(
  tidyverse,
  forecast,
  astsa
)

# Global Tuning Parameters ------------------------------------------------

nlen_ <- 100
v_Mchoice <- 1:10
nsim_ <- 1000
cl_ <- .05

set.seed(1234)


# Simulated Data ----------------------------------------------------------

T_ <- 1000 + nlen_

Simulate_Data <- function(number_obs) {

  dt0_ <- 0
  v_dt <- vector("double", number_obs)

  for (t in 1:number_obs) {

    dt0_ <- .5 + .7 * dt0_ + rnorm(1)
    v_dt[t] <- dt0_

  }

  return(v_dt)
}

v_dt <- Simulate_Data(T_)
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


# Size Computation --------------------------------------------------------

Generate_ARIMA_Simulated <- function(M) {

  v_arima_sim <- arima.sim(m_sim, n = nlen_, innov = rnorm(nlen_, 0, sqrt(a$sigma2)))
  dm_test <- dm.test.bt(v_arima_sim, M=M, cl=cl_)

  return(dm_test)
}

l_arima_sim <- map(1:nsim_, ~ arima.sim(m_sim, n = nlen_, innov = rnorm(nlen_, 0, sqrt(a$sigma2))))
l_dm_test <- map(l_arima_sim, ~ dm.test.bt(., M = 3, cl = cl_))
v_test_statistic <- map_dbl(l_dm_test, ~ pluck(., "stat"))
v_reject <- map_dbl(l_dm_test, ~ pluck(., "rej"))
size_distortion_ <- mean(v_reject) - cl_

paste("Size distortion is", size_distortion_)

# df_sim <- tibble(
#   arima_sim = l_arima_sim,
#   test_statistic = map_dbl(arima_sim, ~ pluck(., "stat")),
#   reject = map_dbl(arima_sim, ~ pluck(., "rej"))
# ) %>%
# print()
#
# size_distortion_ <- mean(df_sim$reject) - cl_




# Delta setup -------------------------------------------------------------

ndel <-  50 #number of deltas
del_tilde <-  10 #largest delta
del_grid <-  seq(from=-del_tilde, to=del_tilde, length.out=ndel)

ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100);
Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than a spectral density at zero

v_oracle_test_statistic <- map_dbl(l_arima_sim, ~ mean(.) / sqrt(Om / nlen_))
v_oracle_p_val <-  2 * stats::pnorm(-abs(v_oracle_test_statistic), mean=0, sd=1)
v_oracle_rej <- v_oracle_p_val < cl_

empirical_size_oracle_ <- mean(v_oracle_rej)

paste("Empirical size of oracle is", empirical_size_oracle_)

c05_star_oracle <- quantile(abs(v_oracle_test_statistic), 1 - cl_)


# size corrected power ----------------------------------------------------

v_delta_sim <- map_d





