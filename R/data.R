# Collection of data

#' Simulated data 1
#'
#' A dataset containing y, f1, f2 based on one of McCracken (2019)'s data generating processes.
#' The name of the specification is "Unconditional-Rolling".
#' R = 175, Rbar = 175, h = 12, P = 75;
#' Note that the DGP is constructed in a way that there is no difference between f1 and f2 in terms of expected quadratic loss.
#' @format A data frame with 75 rows and 3 variables:
#' \describe{
#' \item{y}{Forecast target}
#' \item{f1}{Point forecast based on forecasting model 1}
#' \item{f2}{Point forecast based on forecasting model 2}
#' }
#' @source Simulated data.
"mikedata"
