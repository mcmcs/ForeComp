# ForeComp

This repository contains the **ForeComp** R package for forecasting comparison and equal predictive accuracy testing. The package implements Diebold-Mariano and related procedures, including fixed-smoothing (fixed-b and fixed-m) variants, along with tools for size-power tradeoff analysis.

## Contents

- `R/` package source code implementing forecast evaluation functions
- `data/` packaged datasets used by the package
- `man/` function documentation (`.Rd` files)
- `tests/testthat/` automated tests
- `NEWS.md` release notes and version history

The package code and datasets support forecasting-method comparison and robustness analysis for equal predictive accuracy testing.

## Installation

`ForeComp` requires `R (>= 3.5.0)`.

Install the stable release from CRAN:

```r
install.packages("ForeComp")
```

Install the development version from GitHub:

```r
remotes::install_github("mcmcs/ForeComp")
```

## Main functions

`Plot_Tradeoff(data, f1, f2, y, ...)`

- Computes loss differentials between two forecast series and evaluates how the truncation choice (`M`) affects size distortion and maximum power loss.
- Fits an ARIMA approximation to the loss differential process and uses simulation (`n_sim`) to construct size-power tradeoff results across `m_set`.
- Returns a list where the first element is a `ggplot2` tradeoff figure and the second element is the underlying computed table.

For direct equal-predictive-accuracy testing with specific bandwidth rules:

`ForeComp` supports both fixed-smoothing defaults and baseline alternatives for Bartlett-kernel DM tests:

- `dm.test.bt(d, M = NA, Mopt = ...)` (normal approximation, default `Mopt = 2`)
- `dm.test.bt.fb(d, M = NA, Mopt = ...)` (fixed-b approximation, default `Mopt = 1`)

For both functions, `Mopt` has the same meaning:

- `Mopt = 1` (LLSW): `M = ceiling(1.3 * sqrt(T))`
- `Mopt = 2` (NW 1994): `M = ceiling(4 * (T / 100)^(2/9))`
- `Mopt = 3` (textbook NW / Andrews): `M = ceiling(0.75 * T^(1/3))`
- `Mopt = 4` (CI baseline): `M = floor(sqrt(T))`

where `T = length(d)`.

For EWC fixed-smoothing:

- `dm.test.ewc.fb(d, B = NA, Bopt = 1)` uses `B = floor(0.4 * T^(2/3))` with a lower bound of 1.

## Introductory Paper

Coming soon!
