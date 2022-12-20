
rm(list=ls());

library(forecast);
library(ForeComp)
library(astsa)

# --- Generate fake data
T = 1000+1000;
# T = 100;
dt = matrix(NA, T, 1);
dt0 = 0;
for (t in 1:T){
  dt0 = 0.5 + 0.95*dt0 + rnorm(1);
  dt[t] = dt0;
}
dt = dt[1001:T,,drop=F] # because we are starting from the arbitrary value, we need to discard initial observations

# --- demean
mut= mean(dt);
dt_tilde = dt - mut;

# --- Estimate ARIMA
a = auto.arima(y=dt_tilde, max.p = 12, max.q = 12, stationary=T, ic="aic", seasonal=F, allowmean=F)

# extract info to build a arma model to simulate data
i_ar  = grep("ar", names(a$coef));
i_ma  = grep("ma", names(a$coef));

m_sim = list("ar"=a$coef[i_ar], "ma"=a$coef[i_ma]);

if (length(m_sim$ar)==0){m_sim$ar = 0.0}
if (length(m_sim$ma)==0){m_sim$ma = 0.0}

# --- Size computation
Mchoice = 20;
nlen =20;

nsim = 1000;
# nlen = 50;
mat_rej = matrix(NA, nsim,1);
mat_stat = matrix(NA, nsim,1);
mat_dt = matrix(NA, nlen, nsim); #matrix that stores data (for size calculation)
for (irep in 1:nsim){
  dt_sim = arima.sim(m_sim, n = nlen, innov = rnorm(nlen, 0, sqrt(a$sigma2)));
  rst = dm.test.bt(dt_sim, M = Mchoice, cl = 0.05);
  mat_rej[irep,] = rst$rej;
  mat_stat[irep,] = rst$stat;
  mat_dt[,irep] = dt_sim;

}
print("empirical size");
print(mean(mat_rej));

# find a cut-off (c*) - size-corrected critical value
c05_star = quantile(abs(mat_stat), (1-0.05));

# check empirical rejection rate with c*
nsim = 1000;
mat_rej2 = matrix(NA, nsim,1);
for (irep in 1:nsim){
  dt_sim = mat_dt[,irep,drop=F];
  rst = dm.test.bt(dt_sim, M = Mchoice, cl = 0.05);
  mat_rej2[irep,] = abs(rst$stat) > c05_star; # reject if the statistic is larger than the size-correct crit value in abs
}
print("empirical size");
print(mean(mat_rej2));

# --- Power computation
# H_{0}: E[dt] = 0
del = 0.1;

nsim = 1000;
# nlen = 50;
mat_rej = matrix(NA, nsim,1);  # raw power
mat_rej2 = matrix(NA, nsim,1); # size-adjusted power
mat_stat = matrix(NA, nsim,1);
for (irep in 1:nsim){
  dt_sim = mat_dt[,irep,drop=F];
  dt_sim = dt_sim + del;
  rst = dm.test.bt(dt_sim, M = Mchoice, cl = 0.05);
  mat_rej[irep,] = rst$rej;
  mat_stat[irep,] = rst$stat;

  mat_rej2[irep,] = abs(rst$stat) > c05_star; # test with size-corrected critical value

}
print("empirical rejection rate (power), raw");
print(mean(mat_rej));

# --- Size-adjusted power
print("empirical rejection rate (power), size-corrected");
print(mean(mat_rej2));


# ==============================================================================

# ***to-do
# 1- is one-sided delta range enough?
# 2- loop over M
# 3- make it for all tests (at least three tests we consider in the paper)

# --- Maximum power loss

# some tuning parameters
ndel = 50;
del_tilde = 10;
del_grid = seq(from=-del_tilde, to=del_tilde, length.out = ndel);

cl = 0.05;

# 1: Power of oracle test
# need to find a function that returns long-run variance of y(t) given "a", which is output from auto.arima

# --- Calculate
ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100)
Om = ss$spec[1]; #this is 2*pi*f(0), spectrum at zero rather than spectral density at zero

# --- Oracle test, size
mat_stat_o = matrix(NA,nsim,1);
for (irep in 1:nsim){
  dt_sim = mat_dt[,irep,drop=F];
  dmstat = mean(dt_sim) / sqrt(Om / nlen); #statistic
  mat_stat_o[irep,1] = dmstat;
}

# find a cut-off (c*) - size-corrected critical value
c05_star_o = quantile(abs(mat_stat_o), (1-cl));

# size corrected power
mat_rej_o = matrix(NA,nsim,ndel)
mat_rej2_o = matrix(NA,nsim,ndel)
for (idel in 1:ndel){
  for (irep in 1:nsim){
    dt_sim = mat_dt[,irep,drop=F] + (1/sqrt(nlen)) * sqrt(Om) * del_grid[idel];
    dmstat = mean(dt_sim) / sqrt(Om / nlen); #statistic
    pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation
    rej = pval < cl; #reject decision
    mat_rej_o[irep,idel] = rej;

    mat_rej2_o[irep,idel] = abs(dmstat) > c05_star_o; # test with size-corrected critical value
  }
}
plot(del_grid, apply(mat_rej_o, 2, mean), ylim=range(0,0.5))
lines(del_grid, apply(mat_rej_o, 2, mean)) # raw power
lines(del_grid, apply(mat_rej2_o, 2, mean), col='red') # size -corrected power


# 2: Power given delta
Mchoice = 20;
# --- size
mat_stat = matrix(NA,nsim,1);
for (irep in 1:nsim){
  dt_sim = mat_dt[,irep,drop=F];
  rst = dm.test.bt(dt_sim, M = Mchoice, cl = cl);
  mat_stat[irep,1] = rst$stat;
}
# --- size-corrected crit val
c05_star = quantile(abs(mat_stat), (1-cl));

# --- size-adjusted power
mat_stat= matrix(NA,nsim,ndel)
mat_rej = matrix(NA,nsim,ndel)
mat_rej2 = matrix(NA,nsim,ndel)
for (idel in 1:ndel){


  for (irep in 1:nsim){
    dt_sim = mat_dt[,irep,drop=F] + (1/sqrt(nlen)) * sqrt(Om) * del_grid[idel];
    rst = dm.test.bt(dt_sim, M = Mchoice, cl = cl);
    mat_rej[irep,idel] = rst$rej;
    mat_stat[irep,idel] = rst$stat;
    mat_rej2[irep,idel] = abs(rst$stat) > c05_star; # test with size-corrected critical value
  }

}
plot(del_grid, apply(mat_rej, 2, mean), ylim=range(0,0.4))
lines(del_grid, apply(mat_rej, 2, mean))
lines(del_grid, apply(mat_rej2, 2, mean), col="red")

plot(del_grid, apply(mat_rej2_o, 2, mean), ylim=range(0,1))
lines(del_grid, apply(mat_rej2_o, 2, mean))
lines(del_grid, apply(mat_rej2, 2, mean),col="red")

xxxxxxxxxxx
plot(apply(mat_rej2_o, 2, mean)-apply(mat_rej2, 2, mean))

# ==============================================================================

# this is power spectrum -> how to move to spectral density?

# --- Sanity checks
# --- Comparing the theoretical and simulated long-run variance
nsim   = 1000;
nlen   = 10000;
mat_dt = matrix(NA, nlen, nsim);
for (irep in 1:nsim){
  dt_sim = arima.sim(m_sim, n = nlen, innov = rnorm(nlen, 0, sqrt(a$sigma2)));
  mat_dt[,irep] = dt_sim;
}
ss = arma.spec(ar = m_sim$ar, ma = m_sim$ma, var.noise = a$sigma2, n.freq = 100)

print(ss$spec[1]/nlen)
print(sd(apply(mat_dt, 2, mean))^2)




