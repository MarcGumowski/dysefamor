#' Yield Curves for Canadian Zero-Coupon Bonds
#'
#' A dataset containing daily yields curves for zero-coupon bonds, generated
#' using pricing data for Government of Canada bonds and treasury bills.
#' Each row is a single zero-coupon yield curve. The first column is the date.
#' Each remaining column is a different term to maturity ranging from 0.25 years
#' (column 2) to 30.00 years (column 16). The data are expressed as decimals
#' (e.g. 0.0500 = 5.00\% yield).
#'
#' @format A data frame with 6211 rows and 16 variables: the date and 15
#' different time-to-maturity covariates.

#' @source \url{http://www.bankofcanada.ca/rates/interest-rates/bond-yield-curves/}
#'
"canadianYieldCurves"
