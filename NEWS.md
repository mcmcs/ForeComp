# ForeComp

## ForeComp 1.0.0 (19 February 2026)

- Upgraded `Plot_Tradeoff` with new functionality and stronger robustness checks.
- Exposed additional DM testing procedures as standalone functions for direct use outside `Plot_Tradeoff`.
- Improved Bartlett test handling for large truncation values (`M`) to avoid recycling warnings and stabilize inference.
- Added a full `testthat` suite for core DM-style tests, `loss.diff.p`, input validation, and `Plot_Tradeoff`.
- Updated defaults/behavior checks for bandwidth-selection paths, including local robustness checks for `n_sim = 500`.

## ForeComp 0.9.0 (1 August 2023)

- First public release.
