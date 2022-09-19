
pacman::p_load(
  tidyverse,
  lubridate,
  forecast,
  tictoc
)


Load_Data <- function(dataset) {

   read_csv(paste0(dataset, "_extended.csv")) %>%
    select(X1, starts_with("SPF"), starts_with("NC"), Realiz1) %>%
    mutate(
      year_quarter  = str_replace(X1, ":", "-") %>% yq(.),
      across(-year_quarter, ~ as.double(.))
    ) %>%
    select(-X1) %>%
    select(year_quarter, everything())

}

df_inflation <- Load_Data("PGDP") %>% glimpse()

dt <- df_inflation$Realiz1 - df_inflation$SPFfor_Step1
dt[is.na(dt)] <- 0

# n <- length(dt)
n <- 40

mu_dt <- mean(dt)
dt_tilde <- dt - mu_dt

if(!all.equal(0, mean(dt_tilde))){stop("Mean should be 0")}

arima_model <- auto.arima(y = dt_tilde, max.p = 12, max.q = 12,
                          stationary = TRUE, ic = "aic", seasonal = FALSE, allowmean = FALSE)

Simulate_ARMA <- function(bandwidth) {

  ft_tilde <- arima.sim(model = arima_model$model, n = n)

  wce_b = dm.test.bt.fb(ft_tilde, cl = 0.05, M = bandwidth)$rej; #fixed-b (WCE-B)
  wpe_d = dm.test.wpe.fb(ft_tilde, cl=0.05, M = bandwidth)$rej; # (WPE-D)

  tibble(
    "bandwidth" = bandwidth,
    "wce_b" = wce_b,
    "wpe_d" = wpe_d
  )

}

v_bandwidth <- c(1:10, 20, 50)

set.seed(1234)

# 45 seconds to run
tic()
df_sim <- rerun(1000, map_dfr(v_bandwidth, ~ Simulate_ARMA(.))) %>% # add print statement for which bandwidth we're on
  bind_rows(.) %>%
  group_by(bandwidth) %>%
  summarize(
    across(everything(), ~ mean(.))
  )
toc()

