% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NestRRRini.r
\name{NRRR.ini}
\alias{NRRR.ini}
\title{Generate initial estimators}
\usage{
NRRR.ini(Y, X, r, rx, ry, jx, jy, p, d, n)
}
\arguments{
\item{Y}{response matrix of dimension n-by-jy*d.}

\item{X}{design matrix of dimension n-by-jx*p.}

\item{r}{rank of the local reduced-rank structure.}

\item{rx}{number of latent predictors.}

\item{ry}{number of latent responses.}

\item{jx}{number of basis functions to expand the functional predictor.}

\item{jy}{number of basis functions to expand the functional response.}

\item{p}{number of predictors.}

\item{d}{number of responses.}

\item{n}{sample size.}
}
\value{
The function returns a list:
  \item{Ag}{the estimated U.}
  \item{Bg}{the estimated V.}
}
\description{
This function provides the initial estimators of U and V to initialize the
blockwise coordinate descent algorithm.
}
\examples{
library(NRRR)
simDat <- NRRR.sim(n = 100, ns = 200, nt = 200, r = 5, rx = 3, ry = 3,
                   jx = 15, jy = 15, p = 10, d = 6, s2n = 1, rho_X = 0.5,
                   rho_E = 0, Sigma = "CorrAR")
fit_init <- with(simDat, NRRR.ini(Y = Yest, X = Xest, r = 5,
                 rx = 3, ry = 3, jx = 15, jy = 15,
                 p = 10, d = 6, n = 100))
fit_init$Ag
}
\references{
Liu, X., Ma, S., & Chen, K. (2020).
Multivariate Functional Regression via Nested Reduced-Rank Regularization.
arXiv: Methodology.
}
