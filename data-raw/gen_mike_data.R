# Constructing artificial data
# Mike McCraken (2019)'s DGP, rolling-unconditional

# -----------------------
# House-keeping
rm(list = ls());

workpath = file.path(getwd(), "data-raw");


# -----------------------
# Main settings
R = 175;
Rbar = 175;
h = 12;
P = 75;
seed = 1234;

# -----------------------
# Some Fixed Settings
sig2 = 1;
theta = matrix((0.5)^(0:(h-1)));
T = Rbar + P + h - 1;
sig = sqrt(sig2);

# ---------------------------
# DGP: rolling, unconditional
Gam = matrix(0,nrow=R,ncol=1); #from 0 to (R-1)
for (t in 1:R){
  if (t <= h){
    Gam[t] = sig2 * sum(theta[1:(h-t+1)] * theta[t:h]);
  } else {
    Gam[t] = 0;
  }
}

mu = 1/sqrt(R) * sqrt(Gam[1] + 2*(sum((R-(1:(R-1)))/R * Gam[2:R] )));

# ---------------------------
# function to generate data
genDataEq = function(T,h,mu,theta, sig){
  e0 = matrix(rnorm(n=h,mean=0,sd=sig));
  yt = matrix(0, nrow=T, ncol=1);
  for (t in 1:T){
    yt[t] = mu + sum(e0*theta)
    #e0 = rbind(rnorm(n=1,mean=0,sd=sig), e0[1:(h-1),, drop=F]);
    e0 = rbind(rnorm(n=1,mean=0,sd=sig), e0[(-h), , drop=F]);
  }
  return(yt);
}

# ----------------------------
# Data
yt = genDataEq(T,h,mu,theta,sig);

# ----------------------------
# Evaluation, Rolling

# Rolling prediction and forecast error
mat_y  = matrix(NA, nrow=P, ncol=1);
mat_f1 = matrix(NA, nrow=P, ncol=1);
mat_f2 = matrix(NA, nrow=P, ncol=1);

mat_e1 = matrix(NA, nrow=P, ncol=1);
mat_e2 = matrix(NA, nrow=P, ncol=1);
for (t in (Rbar:(Rbar+P-1))){

  # Prediction (h-step-ahead)
  yhat1 = 0;
  yhat2 = mean(yt[(t-Rbar+1):t]);

  # h-step-ahead
  mat_e1[t-Rbar+1, 1] = yt[t+h] - yhat1;
  mat_e2[t-Rbar+1, 1] = yt[t+h] - yhat2;

  # collect relevant information
  mat_y[(t-Rbar+1), 1] = yt[t+h];
  mat_f1[(t-Rbar+1), 1] = yhat1;
  mat_f2[(t-Rbar+1), 1] = yhat2;

}

# Format it as data.frame
mikedata = data.frame(cbind(mat_y, mat_f1, mat_f2));
names(mikedata) = c("y", "f1", "f2");

# Save data
# usethis::use_data(mikedata, overwrite = TRUE)
usethis::use_data(mikedata)

