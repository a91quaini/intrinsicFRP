#' Test Asset Excess Returns - monthly observations from `07/1963` to `02/2024`
#'
#' Monthly excess returns on the `25` Size/Book-to-Market double sorted portfolios
#' and the `17` industry portfolios from `07/1963` to `02/2024`.
#'
#' @format ## `returns`
#' A data frame with `624` rows and `43` columns:
#' \describe{
#'   \item{Date}{Date in yyyymm format}
#'   \item{...}{Asset excess returns}
#' }
#' @source <https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>
"returns"

#' Factors - monthly observations from `07/1963` to `02/2024`
#'
#' Monthly observations from `07/1963` to `02/2024` of
#' the Fama-French `5` factors and the momentum factor.
#'
#' @format ## `factors`
#' A data frame with `624` rows and `7` columns:
#' \describe{
#'   \item{Date}{Date in yyyymm format}
#'   \item{...}{Factor observations}
#' }
#' @source <https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>
"factors"

#' Risk free - monthly observations from `07/1963` to `02/2024`
#'
#' Monthly observations from `07/1963` to `02/2024` of
#' the US risk free asset.
#'
#' @format ## `risk_free`
#' A data frame with `624` rows and `2` columns:
#' \describe{
#'   \item{Date}{Date in yyyymm format}
#'   \item{...}{risk free observations}
#' }
#' @source <https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>
"risk_free"
