# ForeComp

This repository contains the **ForeComp** R package and replication scripts for a forecasting comparison study. The package implements Diebold-Mariano and related tests for equal predictive accuracy. Datasets and scripts included here reproduce tradeoff plots and robustness tables used in the accompanying paper.

## Contents

- `R/` package source code implementing forecast evaluation functions
- `data/` packaged datasets and `data-raw/` generation inputs/scripts
- Replication scripts under `applications/` (e.g., `applications/application02/ci_replication_minchul.R`, `applications/application03/example_plot_tradeoff.R`)
- Legacy development scripts under `archive/` (kept for record; not used in current dev)

The code and data allow users to replicate the empirical results and explore the performance of different forecasting methods.

## Bandwidth Options

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
