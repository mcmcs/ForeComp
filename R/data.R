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

#' Civilian Unemployment Rate
#'
#' Error Statistics for the Survey of Professional Forecasters for Nominal GNP/GDP
#'
#' @format ## `UNEMP`
#' A data frame with 219 rows and 31 columns:
#' \describe{
#'   \item{country}{Country name}
#'   \item{iso2, iso3}{2 & 3 letter ISO country codes}
#'   \item{year}{Year}
#'   ...
#' }
#' @source <https://www.philadelphiafed.org/-/media/frbp/assets/surveys-and-data/survey-of-professional-forecasters/data-files/unemp/data_spf_error_statistics_unemp_1_aic.xls?la=en&hash=4CAD0B11FEAB6C4D0F30C38965FE3354>
"UNEMP"
