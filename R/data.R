#' Test Asset Excess Returns - monthly observations from `01/1970` to `12/2021`
#'
#' Monthly excess returns on the `25` Size/Book-to-Market double sorted portfolios
#' and the `17` industry portfolios from `01/1970` to `12/2021`.
#'
#' @format ## `returns`
#' A data frame with `624` rows and `43` columns:
#' \describe{
#'   \item{Date}{Date in yyyymm format}
#'   \item{...}{Asset excess returns}
#' }
#' @source <https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>
"returns"

#' Factors - monthly observations from `01/1970` to `12/2021`
#'
#' Monthly observations from `01/1970` to `12/2021` of
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
