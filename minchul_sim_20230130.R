# Minchul tries to replicate

## Setting up

rm(list=ls());
library(forecast)
library(ForeComp)
library(astsa)


# Global tuning parameter
nlen    = 200; # length of time-series data
nsim    = 10000; # number of simulation to compute size and power
cl      = 0.05; # confidence level


# Mchoice = 2;   # bandwidth
M_set = c(seq(from=1, to=10, by=1), seq(from=11, to=100, by=5)); # grid for M
M_set_n = length(M_set); # length of M grid

mat_size_distortion = matrix(NA, M_set_n, 1);
mat_power_loss = matrix(NA, M_set_n, 1);


set.seed(1234); # set seed


# DGP =================================
T = 10000 + nlen;
dt = matrix(NA, T, 1);
dt0 = 0;
for (t in 1:T){
  dt0 = 0.0 + 0.7*dt0 + rnorm(1);
  dt[t] = dt0 + 0*rnorm(1);
}
dt = dt[1001:T, , drop=F]; # because we are starting from the arbitrary value ...
# we need to discard initial observationas

mut = mean(dt); # demean
dt_tilde = dt - mut;


# --- Estimate ARIMA and extract information
a = auto.arima(y=dt_tilde, max.p = 12, max.q=12, stationary=T, ic="aic", seasonal=F, allowmean=F);
i_ar = grep("ar", names(a$coef));
i_ma = grep("ma", names(a$coef));
m_sim = list("ar"=a$coef[i_ar], "ma"=a$coef[i_ma]);

if (length(m_sim$ar)==0){m_sim$ar=0.0}
if (length(m_sim$ma)==0){m_sim$ma=0.0}

# --- Generate data for size and power computation (this loop alos pre-computes long-run variance)
Mmax = min(c(nlen-1,200)); # maximum possible M considered in this experimen
mat_dt   = matrix(NA, nlen, nsim); # matrix that stores data (for size calculation)
mat_dtm  = matrix(NA, 1, nsim);    # matrix that stores mean of data
mat_acf  = matrix(NA, (Mmax+1), nsim); # matrix that stores acf
for (irep in 1:nsim){

  # data
  dt_sim = arima.sim(m_sim, n=nlen, innov = rnorm(nlen, 0, sqrt(a$sigma2)));
  mat_dt[,irep] = dt_sim;
  mat_dtm[,irep] = mean(dt_sim);
  # autocovariance matrix
  d.cov = stats::acf(dt_sim, lag.max = (Mmax), type="covariance", plot=FALSE, demean=TRUE)$acf[,,1]; #should I compute ACF under the null regardless? Then, this should be deman = FALSE
  mat_acf[,irep] = d.cov;
}


# --- Setting for trade-off figure
# ndel = 50; # number of deltas
# del_tilde = 10; # largest delta
# del_grid = seq(from=-del_tilde, to=del_tilde, length.out=ndel);
del_grid = seq(from=0, to=10, by=0.25);
ndel = length(del_grid);



# --- Loop over M set
for (iM in 1:M_set_n){

  Mchoice = M_set[iM]; #our choice of M for this iteration


  # --- Oracle ---

  # --- Power for oracle

  # long-run variance
  ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100);
  Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than a spectral density at zero

  # # Oracle test
  # mat_stat_o = matrix(NA, nsim, 1);
  # mat_rej_o = matrix(NA, nsim, 1);
  # for (irep in 1:nsim){
  #   dmstat = mean(mat_dtm[,irep]) / sqrt(Om / nlen); #Oracle's statistic
  #   pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
  #   rej = pval < cl; #reject decision
  #   mat_stat_o[irep, 1] = dmstat;
  #   mat_rej_o[irep,1] = rej;
  # }
  # # print("empirical size of oracle");
  # # print(mean(mat_rej_o));
  #
  #
  # # --- find a cut-off (c*) to get size-corrected critical value
  # c05_star_o = quantile(abs(mat_stat_o), (1-cl));
  #
  # # size corrected power
  # mat_rej_o  = matrix(NA, nsim, ndel);
  # mat_rej2_o = matrix(NA, nsim, ndel);
  # for (idel in 1:ndel){
  #   for (irep in 1:nsim){
  #     dmstat = (mat_dtm[,irep] + (1/sqrt(nlen)) *sqrt(Om)*del_grid[idel]) / sqrt(Om / nlen); #statistic
  #     pval = 2*stats::pnorm(-abs(dmstat), mean=0, sd=1); # p-val from normal approximation
  #     rej = pval < cl; # reject decision
  #     mat_rej_o[irep, idel] = rej; # this is for raw power
  #     mat_rej2_o[irep, idel] = abs(dmstat) > c05_star_o; # for size-corrected power
  #   }
  # }



  # --- Standard test


  # --- Size computation for DM-NW
  M = Mchoice;
  mat_stat = matrix(NA, nsim, 1);
  mat_rej  = matrix(NA, nsim, 1);
  mat_d.var = matrix(NA, nsim, 1);
  for (irep in 1:nsim){
    d.cov  = mat_acf[1:(M+1),irep];
    d.var  = ( d.cov[1] + 2*sum( (1 - ((1:M)/M) ) * d.cov[-1] ) ) / nlen;

    # dt_sim = mat_dt[,irep];
    # d.cov = stats::acf(dt_sim, lag.max = (M), type="covariance", plot=FALSE, demean=FALSE)$acf[,,1]; #should I compute ACF under the null regardless? Then, this should be deman = FALSE
    # d.var  = ( d.cov[1] + 2*sum( (1 - ((1:M)/M) ) * d.cov[-1] ) ) / nlen;

    dmstat = mean(mat_dtm[,irep]) / sqrt(d.var);
    pval   = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
    rej    = pval < cl; #reject decision

    # store results
    mat_stat[irep,] = dmstat;
    mat_rej[irep, ] = rej;
    mat_d.var[irep,] = d.var;
  }
  size_distortion = mean(mat_rej) - cl;


  # --- Power of a standard test (DM-NW)
  # --- size-corrected crit val
  c05_star = quantile(abs(mat_stat), (1-cl));

  # --- size-corrected power
  mat_stat = matrix(NA, nsim, ndel);
  mat_rej  = matrix(NA, nsim, ndel);
  mat_rej2 = matrix(NA, nsim, ndel);
  for (irep in 1:nsim){

    m_i = mat_dtm[,irep];
    s_i = sqrt(mat_d.var[irep,]);

    for (idel in 1:ndel){

      # DM test with NW
      dmstat = (m_i + (1/sqrt(nlen)) *sqrt(Om)*del_grid[idel]) / s_i; #statistic
      pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
      rej = pval < cl; #reject decision

      mat_rej[irep, idel] = rej; #for raw power
      mat_stat[irep, idel] = dmstat;
      mat_rej2[irep, idel] = abs(dmstat) > c05_star; # for size-corrected power
    }
  }



  #======================================

  # --- another Oracle ---
  grid = seq(from=0, to=5, by=0.05);
  samp = del_grid;
  pow_gau = pnorm(-1.96+grid) + pnorm(-1.96-grid);
  powinterp = approx(samp, apply(mat_rej2, 2, mean), grid)$y; #interpolated power
  # plot(grid, approx(samp, apply(mat_rej2, 2, mean), grid)$y)
  # lines(grid, pow_gau)
  max_power_loss = max(pow_gau-powinterp)


  # --- Maximum power loss
  # max_power_loss = max( apply(mat_rej2_o, 2, mean) - apply(mat_rej2, 2, mean) );

  # --- Collect results
  print(paste0("M = ", iM, " / ", M_set_n));
  mat_size_distortion[iM] = size_distortion;
  mat_power_loss[iM] = max_power_loss;

} #end of iM iteration

# final figure
plot(mat_size_distortion, mat_power_loss)

# # --- plotting (two powers: raw versus size-adjusted power)
# plot(del_grid, apply(mat_rej, 2, mean), ylim=range(0, 1.0))
# lines(del_grid, apply(mat_rej2, 2, mean))
#
#
# #  --- Power curves
# plot(del_grid, apply(mat_rej2_o, 2, mean), ylim=range(0,1));
# lines(del_grid, apply(mat_rej2_o, 2, mean));
# lines(del_grid, apply(mat_rej2, 2, mean), col = "red");






