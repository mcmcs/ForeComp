# Functions to perform tests
# 9/15/2019

## To build a documentation (in help window)
## devtools::document()

## To build a manual (in pdf)
## devtools::build_manual(path=getwd())

# ----------------------------------------------------------------------
# Two-sided tests
# ----------------------------------------------------------------------


#' Diebold-Mariano Test (with an original recommendation)
#'
#' This function performs the Diebold-Mariano test with their original recommendation.
#' Let \eqn{d_{t}}{d(t)} be a sequence of loss differential, \eqn{t=1,2,...,T}.
#' Then, the function performs a statistical test for the following null hypothesis:
#' \deqn{E[d_{t}] = 0}{E[d(t)] = 0}
#' This function follows the original recommendation by Diebold and Mariano (1995),
#' where the long-run variance is estimated using the rectangular kernel truncated at
#' \eqn{(h-1)}. This function performs a two-sided test. Normal approximation is used to calculate critical values.
#' @param d loss differential
#' @param h h-step-ahead forecast (default = 1, i.e., comparing one-step-ahead forecasts)
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
#' @keywords internal
#' @noRd


dm.test.r = function(d, h=1, cl = 0.05){

  n = length(d);
  d.cov = stats::acf(d, lag.max = (h-1), type="covariance", plot=FALSE)$acf[,,1];
  d.var = sum(c(d.cov[1], 2*d.cov[-1]))/n;


  dmstat = mean(d)/sqrt(d.var);

  pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val

  rej = pval < cl; #(cl)-level testing, rej-decision

  if (d.var<0){
    # print("dm.eq: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}


#' Diebold-Mariano Test (Modified-DM)
#'
#' Diebold-Mariano Test (Modified-DM). Finite-sample modification to the original DM's test.
#' This is a two-sided test.
#'
#' This function is based on \link[forecast]{dm.test} in the "forecast" package on CRAN.
#'
#' @param d loss differential
#' @param h h-step-ahead forecast (default = 1, i.e., one-step-ahead forecasts)
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.r.m = function(d, h=1, cl = 0.05){

  n = length(d);
  d.cov = stats::acf(d, lag.max = (h-1), type="covariance", plot=FALSE)$acf[,,1];
  d.var = sum(c(d.cov[1], 2*d.cov[-1]))/n;


  dmstat = mean(d)/sqrt(d.var);

  # Finite-sample modification
  k = ((n + 1 - 2*h + (h/n) * (h - 1))/n)^(1/2);
  dmstat = dmstat * k; #finite sample correction

  # Another finite-sample modification (to use t-dist)
  pval = 2 * stats::pt(-abs(dmstat), df = (n-1) ); #p-val

  rej = pval < cl; #(cl)-level testing, rej-decision

  if (d.var<0){
    # print("dm.eq: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}


# ---------------------------------------------------------------
#' Diebold-Mariano Test (Bartlett kernel, normal approximation)
#'
#' Diebold-Mariano Test (Bartlett kernel, normal approximation). This is a two-sided test.
#'
#' @param d loss differential
#' @param M truncation parameter for the Bartlett kernel (if \code{M = NA}, then \code{Mopt = 2} by default)
#' @param Mopt option for optimal bandwidth, 1 if Lazarus et al. (2018), 2 if Newey and West (1994)
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.bt = function(d, M = NA, Mopt = NA, cl = 0.05){

  n = length(d);

  # Bandwidth selection (if M is not specified)
  # For Bartlett Kernel with normal approximation, NW (1994) is a default
  if (is.na(M)){
    if (is.na(Mopt)){
      Mopt = 2; #Newey and West (1994)
    }

    if (Mopt == 1){
      M = floor(1.3*sqrt(n)) + 1; #Lazarus et al. (2018)
    } else if (Mopt == 2) {
      M = floor(4*(n/100)^(2/9)) + 1; #Newey and West (1994);
    }
  }

  # Long-run variance estimation
  d.cov = stats::acf(d, lag.max = M, type="covariance",plot=FALSE)$acf[,,1];
  d.var = ( d.cov[1] + 2* sum( ( 1 - ( (1:M)/M ) ) * d.cov[-1] ) ) / n;

  dmstat = mean(d) / sqrt(d.var);

  pval = 2 * stats::pnorm(-abs(dmstat), mean=0, sd=1); #p-val based on normal approximation

  rej = pval < cl; #reject decision

  if (d.var<0){
    # Note: in principle, it cannot be negative with prob. one.
    # print("dm.eq.bt: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}

# ---------------------------------------------------------------
#' Diebold-Mariano Test (Bartlett kernel, fixed-b approximation)
#'
#' Diebold-Mariano Test (Bartlett kernel, fixed-b approximation). This is a two-sided test.
#'
#' @param d loss differential
#' @param M truncation parameter for the Bartlett kernel (if \code{M = NA}, then \code{Mopt = 1} by default)
#' @param Mopt option for optimal bandwidth, 1 if Lazarus et al. (2018), 2 if Newey and West (1994)
#' @param cl confidence level (default = 0.05, i.e., 5\%), Only 0.05 (5\%) or 0.10 (10\%) are allowed.
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.bt.fb = function(d, M = NA, Mopt = NA, cl = 0.05){

  n = length(d);

  # Bandwidth selection (if M is not specified)
  # For Bartlett Kernel with fixed-b approximation, Lazarus et al. is a default
  if (is.na(M)){
    if (is.na(Mopt)){
      Mopt = 1; #Lazarus et al. (2018)
    }

    if (Mopt == 1){
      M = floor(1.3*sqrt(n)) + 1; #Lazarus et al. (2018)
    } else if (Mopt == 2) {
      M = floor(4*(n/100)^(2/9)) + 1; #Newey and West (1994);
    }
  }

  # Long-run variance estimation
  d.cov = stats::acf(d, lag.max = M, type="covariance",plot=FALSE)$acf[,,1];
  d.var = ( d.cov[1] + 2* sum( ( 1 - ( (1:M)/M ) ) * d.cov[-1] ) ) / n;

  dmstat = mean(d) / sqrt(d.var);

  # Kiefer and Vogelsang (2005)
  b = M/n;
  #crit950 = 1.6449 + 2.1859*b + 0.3142*b^2 -0.3427*b^3; #0.950 quantile
  #crit975 = 1.9600 + 2.9694*b + 0.4160*b^2 -0.5324*b^3; #0.975 quantile
  if (cl == 0.05){
    crit = 1.9600 + 2.9694*b + 0.4160*b^2 -0.5324*b^3; #0.975 quantile
  } else if (cl == 0.10){
    crit = 1.6449 + 2.1859*b + 0.3142*b^2 -0.3427*b^3; #0.950 quantile
  } else {
    print("dm.eq.bt.fb: cl has to be 0.05 or 0.10 ... function exited without testing");
    outls = list();
    return(outls);
  }

  # rejection decision
  rej = abs(dmstat) > crit; #reject decision

  if (d.var<0){
    # Note: in principle, it cannot be negative with prob. one.
    # print("dm.eq.bt.fb: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;

  return(outls);
}


# ---------------------------------------------------------------
#' Diebold-Mariano Test (EWC, fixed-b approximation)
#'
#' Diebold-Mariano Test (EWC, fixed-b approximation). This is a two-sided test.
#'
#' @param d loss differential
#' @param B truncation parameter for the EWC long-run variance estimator (if \code{B = NA}, then \code{Bopt = 1} by default)
#' @param Bopt option for optimal bandwidth, 1 if Lazarus et al. (2018)' recommendation
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.ewc.fb = function(d, B = NA, Bopt = NA, cl = 0.05){

  n = length(d);

  # Bandwidth selection (if M is not specified)
  # For Bartlett Kernel with fixed-b approximation, Lazarus et al. is a default
  if (is.na(B)){
    if (is.na(Bopt)){
      Bopt = 1; #Lazarus et al. (2018)
    }

    if (Bopt == 1){
      B = floor(0.4*n^(2/3)); #Lazarus et al. (2018)
    }
  }

  # EWC estimator
  j_ind = matrix(1:B);

  t_ind = ( t(matrix(1:n)) - 1/2 ) / n;

  Lhat = sqrt(2/n) * ( cos(pi*j_ind %*% t_ind) %*% d ); #eqn (10) of Lazarus et al. (2018)
  Ohat = crossprod(Lhat) / B; #eqn (11)

  d.var = Ohat / n;
  dmstat = mean(d) / sqrt(d.var);

  pval = 2 * stats::pt(-abs(dmstat), df = B); #p-val


  rej = pval < cl; #reject decision

  if (d.var<0){
    # Note: in principle, it cannot be negative with prob. one.
    # print("dm.eq.ewc: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}


# ---------------------------------------------------------------
#' Diebold-Mariano Test (WPE, fixed-m approximation)
#'
#' Diebold-Mariano Test (WPE, fixed-m approximation). This is a two-sided test. See Coroneo and Iacone (2020)
#'
#' @param d loss differential
#' @param M truncation parameter for the WPE long-run variance estimator (if \code{M = NA}, then \code{Mopt = 1} by default)
#' @param Mopt option for optimal bandwidth, 1 if Coroneo and Iacone's default value (M = floor(T^(1/3)))
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.wpe.fb = function(d, M = NA, Mopt = NA, cl = 0.05){

  n = length(d);

  # Bandwidth selection (if M is not specified)
  if (is.na(M)){
    if (is.na(Mopt)){
      Mopt = 1; # Coroneo and Iacone's default
    }

    if (Mopt == 1){
      M = floor(n^(1/3)); #Lazarus et al. (2018)
    }
  }

  # WPE as in Coroneo and Iacone (2020)
  lamj = matrix(2*pi/n * (1:M));
  gam2 = 2*pi*mean( abs( 1/(sqrt(2*pi*n)) * ( exp( -1i * ( lamj %*% t( matrix(1:n)) ) ) %*% d ) )^2 ); # Eqn(8)
  dmstat = sqrt(n) * mean(d)/sqrt(gam2) #Eqn(9)
  pval = 2 * stats::pt(-abs(dmstat), df = 2*M); #p-val
  rej = pval < cl; #reject decision

  if (gam2<0){
    # Note: in principle, it cannot be negative with prob. one.
    # print("dm.eq.ewc: negative variance ... ");
    rej = NA; #negative variance
  }

  #output
  outls = list();
  outls$rej = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}

# ---------------------------------------------------------------
#' Ibragimov and Muller (2010)
#'
#' t-Statistic based HAR-inference by Ibragimov and Muller (2010).
#'
#' @param d loss differential
#' @param q number of blocks
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.im = function(d, q = 2, cl = 0.05){

  n = length(d);

  # b = n/q; #block size (q should be multiple of n, it does not have to be, but for an ease of computation ... )

  # poorman's way of getting out from the case in which n is not a multiple of q
  b = ceiling(n/q);
  D = rep(NA,(b*q));
  D[1:n] = d;

  m = apply(matrix(D, nrow = b, ncol = q), 2, mean, na.rm=T);

  # eqn (17) in CRS
  mbar   = mean(m);
  s2bar  = sum((m-mbar)^2)/(q-1);
  dmstat = mbar / sqrt(s2bar/q);

  # critical value from t-distribution
  pval = 2 * stats::pt(-abs(dmstat), df = (q-1)); #p-val
  rej = pval < cl; #reject decision

  #output
  outls = list();
  outls$rej  = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}


# ---------------------------------------------------------------
#' CNR (2017), |t|-test statitic
#'
#' Randomization test based on asymptotic symmetry, |t|-test statistic
#'
#' @param d loss differential
#' @param q number of blocks
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @param R number of randomized test stats
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.cnr.t = function(d, q = 2, cl = 0.05){

  n = length(d);

  # b = n/q; #block size (q should be multiple of n, it does not have to be, but for an ease of computation ... )

  # poorman's way of getting out from the case in which n is not a multiple of q
  b = ceiling(n/q);
  D = rep(NA,(b*q));
  D[1:n] = d;

  # ---
  # t-version, eqn (17) in CRS
  Sj     = sqrt(n) * as.matrix(apply(matrix(D, nrow = b, ncol = q), 2, mean, na.rm=T));
  Sbar   = mean(Sj);
  Sigbar = sum( (Sj-Sbar)^2 ) / (q-1);
  dmstat = abs(Sbar) / sqrt(Sigbar / q);
  # sign transformation of t-stats
  M = 2^q;
  G = as.matrix(expand.grid(Map(c, rep(-1,q), rep(1,q))));
  Sj_k = G * rep(Sj, rep.int(M, q));
  Sbar_k = apply(Sj_k, 1, mean);
  Sigbar_k = apply( ( Sj_k - Sbar_k %*% matrix(1, nrow = 1, ncol = q) )^2, 1, sum) / (q-1);
  dmstat_k = abs(Sbar_k) / sqrt(Sigbar_k / q);
  # ---

  # # ---
  # # wald-version, eqn (16) in CRS
  # Sj = sqrt(n) * as.matrix(apply(matrix(D, nrow = b, ncol = q), 2, mean, na.rm=T));
  # Sbar = mean(Sj);
  # Sigbar = as.double(t(Sj)%*%Sj / q);
  # dmstat = q * ( (Sbar)^2 ) / Sigbar; #Wald-type test stat
  # # sign transformed of wald-stats
  # # note in d=1, there is no effect of sign change on Sigbar
  # M = 2^q;
  # G = as.matrix(expand.grid(Map(c, rep(-1,q), rep(1,q))));
  # Sbar_k   = G%*%Sj / q;
  # Sigbar_k = Sigbar; # note in d=1, there is no effect of sign change on Sigbar
  # dmstat_k = q*(Sbar_k^2) / Sigbar_k;
  # # ---

  # ---
  # decision-rule (1 if reject)
  dmstat_k = sort(dmstat_k);
  k = ceiling(M*(1-cl));
  if (dmstat == dmstat_k[k]){

    Mp = sum(dmstat_k > dmstat_k[k]);
    M0 = sum(dmstat_k == dmstat_k[k])

    a = (M*cl - Mp) / M0;
    rej = (1 == rbinom(n=1, size=1, prob=a)); #randomization

  } else {
    rej = (dmstat > dmstat_k[k]);
  }

  # p-value
  pval = mean(dmstat >= dmstat_k); #eqn (9) in CNR

  #output
  outls = list();
  outls$rej  = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}

# ---------------------------------------------------------------
#' CNR (2017), Wald-statistic
#'
#' Randomization test based on asymptotic symmetry, Wald-statistic
#'
#' @param d loss differential
#' @param q number of blocks
#' @param cl confidence level (default = 0.05, i.e., 5\%)
#' @param R number of randomized test stats
#' @return This function returns a class with the following elements
#' \itemize{
#' \item \code{rej} is a T/F value. TRUE (reject), FALSE (accept)
#' \item \code{stat} is a test statistic
#' \item \code{pval} is an associated p-value
#' }
#' @author Minchul Shin
#' @export
dm.test.cnr.w = function(d, q = 2, cl = 0.05){

  n = length(d);

  # b = n/q; #block size (q should be multiple of n, it does not have to be, but for an ease of computation ... )

  # poorman's way of getting out from the case in which n is not a multiple of q
  b = ceiling(n/q);
  D = rep(NA,(b*q));
  D[1:n] = d;

  # # ---
  # # t-version, eqn (17) in CRS
  # Sj     = sqrt(n) * as.matrix(apply(matrix(D, nrow = b, ncol = q), 2, mean, na.rm=T));
  # Sbar   = mean(Sj);
  # Sigbar = sum( (Sj-Sbar)^2 ) / (q-1);
  # dmstat = abs(Sbar) / sqrt(Sigbar / q);
  # # sign transformation of t-stats
  # M = 2^q;
  # G = as.matrix(expand.grid(Map(c, rep(-1,q), rep(1,q))));
  # Sj_k = G * rep(Sj, rep.int(M, q));
  # Sbar_k = apply(Sj_k, 1, mean);
  # Sigbar_k = apply( ( Sj_k - Sbar_k %*% matrix(1, nrow = 1, ncol = q) )^2, 1, sum) / (q-1);
  # dmstat_k = abs(Sbar_k) / sqrt(Sigbar_k / q);
  # # ---

  # ---
  # wald-version, eqn (16) in CRS
  Sj = sqrt(n) * as.matrix(apply(matrix(D, nrow = b, ncol = q), 2, mean, na.rm=T));
  Sbar = mean(Sj);
  Sigbar = as.double(t(Sj)%*%Sj / q);
  dmstat = q * ( (Sbar)^2 ) / Sigbar; #Wald-type test stat
  # sign transformed of wald-stats
  # note in d=1, there is no effect of sign change on Sigbar
  M = 2^q;
  G = as.matrix(expand.grid(Map(c, rep(-1,q), rep(1,q))));
  Sbar_k   = G%*%Sj / q;
  Sigbar_k = Sigbar; # note in d=1, there is no effect of sign change on Sigbar
  dmstat_k = q*(Sbar_k^2) / Sigbar_k;
  # ---

  # ---
  # decision-rule (1 if reject)
  dmstat_k = sort(dmstat_k);
  k = ceiling(M*(1-cl));
  if (dmstat == dmstat_k[k]){

    Mp = sum(dmstat_k > dmstat_k[k]);
    M0 = sum(dmstat_k == dmstat_k[k])

    a = (M*cl - Mp) / M0;
    rej = (1 == rbinom(n=1, size=1, prob=a)); #randomization

  } else {
    rej = (dmstat > dmstat_k[k]);
  }

  # p-value
  pval = mean(dmstat >= dmstat_k); #eqn (9) in CNR

  #output
  outls = list();
  outls$rej  = rej;
  outls$stat = dmstat;
  outls$pval = pval;

  return(outls);
}


# # eqn (17) in CRS
# mbar   = mean(m);
# s2bar  = sum((m-mbar)^2)/(q-1);
# dmstat = mbar / sqrt(s2bar/q) ;
#
# # randomization
# # rr = 1000;
# p_dmstat = rep(NA,times=R);
#
# for (rrind in 1:R) {
#   p_m = m*sample(c(-1,1), size=q, replace=TRUE);
#   p_mbar   = mean(p_m);
#   p_s2bar  = sum((p_m-p_mbar)^2)/(q-1);
#   p_dmstat[rrind] = p_mbar / sqrt(p_s2bar/q);
# }
#
# # pval = 2*mean(p_dmstat < -abs(dmstat));
# # pval = mean( (p_dmstat < -abs(dmstat))&(p_dmstat > abs(dmstat)) ) ;
# pval = mean( abs(p_dmstat) <= abs(dmstat) );
#
#
# rej = pval < cl; #reject decision
