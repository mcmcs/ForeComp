
pacman::p_load(tidyverse)

library(ForeComp)

# load data


# ==========================================
# example (1-step ahead (0-horizon in the paper) ); spf versus nc), Full sample

gdp <- read_csv("RGDP_US.csv", col_types = "cddddddddddd")

# gdp = na.omit(gdp); # I don't know why we have a missing in the realization, but let us omit it to see if this replicates the table

e1 = gdp$Realiz1 - gdp$NCfor_Step1; #forecast error from forecaster 1
e2 = gdp$Realiz1 - gdp$SPFfor_Step1; #forecast error from forecaster 2
d  = e1^2 - e2^2; #squared loss differential
d  = na.omit(d);
n = length(d);

wce_dm = dm.test.r(d, cl = 0.05, h=1) #dm original
wce_b = dm.test.bt.fb(d, cl = 0.05, M = floor(sqrt(n))); #fixed-b
wce_b_opt = dm.test.bt.fb(d, cl = 0.05, Mopt=1); #fixed-b with the recent recommendation

wce_b = dm.test.bt.fb(d, cl = 0.1, M = floor(sqrt(n))); #fixed-b
wce_b_opt = dm.test.bt.fb(d, cl = 0.1, Mopt=1); #fixed-b with the recent recommendation


# ==========================================
# example (1-step ahead (0-horizon in the paper) ); spf versus nc)
# Sample range: 2006:Q1-2016Q4
# -> This is range that does not contain missing values -> we get the same number

gdp <- read_csv("RGDP_US.csv", col_types = "cddddddddddd")

e1 = gdp$Realiz1 - gdp$NCfor_Step1; #forecast error from forecaster 1
e2 = gdp$Realiz1 - gdp$SPFfor_Step1; #forecast error from forecaster 2
d  = e1^2 - e2^2; #squared loss differential
d  = d[81:120]
n = length(d);

# 5%
wce_dm = dm.test.r(d, cl = 0.05, h=1) #dm original (WCE-DM)
wce_b = dm.test.bt.fb(d, cl = 0.05, M = floor(sqrt(n))); #fixed-b (WCE-B)
wpe_d = dm.test.wpe.fb(d, cl=0.05, M = floor(n^(1/3))); # (WPE-D)


print("========================================")
print("First column for 2006Q1-2016Q4")
print(matrix(c(wce_dm$stat, wce_b$stat, wce_b_opt$stat)))
print("========================================")


# wce_b_opt = dm.test.bt.fb(d, cl = 0.05, Mopt=1); #fixed-b with the modern recommendation value

## 10%
# wce_b = dm.test.bt.fb(d, cl = 0.1, M = floor(sqrt(n))); #fixed-b
# wce_b_opt = dm.test.bt.fb(d, cl = 0.1, Mopt=1); #fixed-b with the recent recommendation
# wce_d_opt = dm.test.ewc.fb(d, cl=0.1, Bopt=1);
# wce_d_opt = dm.test.ewc.fb(d, cl=0.1, B = floor(n^(1/3)));





