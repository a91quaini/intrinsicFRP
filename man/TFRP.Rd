% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfrp.R
\name{TFRP}
\alias{TFRP}
\title{Compute tradable factor risk premia}
\usage{
TFRP(returns, factors, include_standard_errors = FALSE, check_arguments = TRUE)
}
\arguments{
\item{returns}{\verb{n_observations x n_returns}-dimensional matrix of test asset
excess returns.}

\item{factors}{\verb{n_observations x n_factors}-dimensional matrix of factors.}

\item{include_standard_errors}{boolean \code{TRUE} if you want to compute the
tradable factor risk premia HAC standard errors; \code{FALSE} otherwise. Default
is \code{FALSE}.}

\item{check_arguments}{boolean \code{TRUE} if you want to check function arguments;
\code{FALSE} otherwise. Default is \code{TRUE}.}
}
\value{
a list containing \code{n_factors}-dimensional vector of tradable factor
risk premia in \code{"risk_premia"}; if \code{include_standard_errors=TRUE}, then
it further includes \code{n_factors}-dimensional vector of tradable factor risk
premia standard errors in \code{"standard_errors"}.
}
\description{
Computes tradable factor risk premia from data on factors and
test asset excess returns. Optionally computes the corresponding
heteroskedasticity and autocorrelation robust standard errors using the
Newey-West (1994) plug-in procedure to select the number of relevant lags,
i.e., \code{n_lags = 4 * (n_observations/100)^(2/9)}.
}
\examples{
# import package data on 15 risk factors and 42 test asset excess returns
factors = intrinsicFRP::factors[,-1]
returns = intrinsicFRP::returns[,-1]

# compute tradable factor risk premia and their standard errors
tfrp = TFRP(returns, factors, include_standard_errors = TRUE)

}