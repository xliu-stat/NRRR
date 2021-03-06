% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NestRRR.select.r
\name{NRRR.ic}
\alias{NRRR.ic}
\title{Select ranks with information criterion}
\usage{
NRRR.ic(Y, X, Ag0 = NULL, Bg0 = NULL,
        jx, jy, p, d, n, maxiter = 300, conv = 1e-4,
        method = c('RRR','RRS')[1],
        lambda = 0, ic = c("BIC","BICP","AIC","GCV")[1],
        dimred = c(TRUE,TRUE,TRUE), rankfix = NULL,
        xrankfix = NULL, yrankfix = NULL,
        lang = c('R', 'Rcpp')[1])
}
\arguments{
\item{Y}{response matrix of dimension n-by-jy*d.}

\item{X}{design matrix of dimension n-by-jx*p.}

\item{Ag0}{an initial estimator of matrix U. If \code{Ag0 = NULL} then generate it
by \code{\link{NRRR.ini}}. Default is NULL. If \code{lang = 'Rcpp'},
then \code{Ag0} is automatically generated by \code{\link{NRRR.ini}}.}

\item{Bg0}{an initial estimator of matrix V, if \code{Bg0 = NULL} then generate it
by \code{\link{NRRR.ini}}. Default is NULL. If \code{lang = 'Rcpp'},
then \code{Bg0} is automatically generated by \code{\link{NRRR.ini}}.}

\item{jx}{number of basis functions to expand the functional predictor.}

\item{jy}{number of basis functions to expand the functional response.}

\item{p}{number of predictors.}

\item{d}{number of responses.}

\item{n}{sample size.}

\item{maxiter}{the maximum iteration number of the
blockwise coordinate descent algorithm. Default is 300.}

\item{conv}{the tolerance level used to control the convergence of the
blockwise coordinate descent algorithm. Default is 1e-4.}

\item{method}{'RRR' (default): no additional ridge penalty; 'RRS': add an
additional ridge penalty.}

\item{lambda}{the tuning parameter to control the amount of ridge
penalization. It is only used when \code{method = 'RRS'}.
Default is 0.}

\item{ic}{the user-specified information criterion. Four options are available,
i.e., BIC, BICP, AIC, GCV.}

\item{dimred}{a vector of logical values to decide whether to use the selected information criterion
do rank selection on certain dimensions. TRUE means the rank is selected
by the selected information criterion. If \code{dimred[1] = FALSE}, r is
provided by \code{rankfix} or \eqn{min(jy*d, rank(X))};
If \code{dimred[2] = FALSE}, rx equals to \code{xrankfix} or p; If \code{dimred[3] = FALSE},
ry equals to \code{yrankfix} or d. Default is \code{c(TRUE, TRUE, TRUE)}.}

\item{rankfix}{a user-provided value of r when \code{dimred[1] = FALSE}. Default is NULL
which leads to \eqn{r = min(jy*d, rank(X))}.}

\item{xrankfix}{a user-provided value of rx when \code{dimred[2] = FALSE}. Default is NULL
which leads to \code{rx = p}.}

\item{yrankfix}{a user-provided value of ry when \code{dimred[3] = FALSE}. Default is NULL
which leads to \code{ry = d}.}

\item{lang}{'R' (default): the R version function is used; 'Rcpp': the Rcpp
version function is used.}
}
\value{
The function returns a list:
  \item{Ag}{the estimated U.}
  \item{Bg}{the estimated V.}
  \item{Al}{the estimated A.}
  \item{Bl}{the estimated B.}
  \item{C}{the estimated coefficient matrix C.}
  \item{df}{the estimated degrees of freedom of the NRRR model.}
  \item{sse}{the sum of squared errors of the selected model.}
  \item{ic}{a vector containing values of BIC, BICP, AIC, GCV of the selected model.}
  \item{rank}{the estimated r.}
  \item{rx}{the estimated rx.}
  \item{ry}{the estimated ry.}
}
\description{
This function selects the optimal ranks \code{(r, rx, ry)} using a user-specified
information criterion. The degrees of freedom is estimated by the number of
free parameters in the model. The blockwise coordinate descent algorithm is used to fit
the model with any combinations of \code{(r, rx, ry)}.
}
\details{
Denote \eqn{\hat C(r, rx, ry)} as the estimator of the coefficient matrix
with the rank values fixed at some \eqn{(r, rx, ry)} and write the sum of
squared errors as \eqn{SSE(r, rx, ry)=||Y - X\hat C(r, rx, ry)||_F^2}. We
define
\deqn{BIC(r, rx, ry) = n*d*jy*log(SSE(r, rx, ry)/(n*d*jy)) + log(n*d*jy)*df(r, rx, ry),}
where \eqn{df(r, rx, ry)} is the effective degrees of freedom and is estimated
by the number of free parameters
\deqn{\hat df(r, rx, ry) = rx*(rank(X)/jx - rx) + ry*(d - ry) + (jy*ry + jx*rx - r)*r.}
We can define other information criteria in a similar way. With the defined BIC,
a three-dimensional grid search procedure of the rank
values is performed, and the best model is chosen as the one with the
smallest BIC value. Instead of a nested rank selection method, we apply a
one-at-a-time selection approach. We first set \eqn{rx = p, ry = d}, and
select the best local rank \eqn{\hat r} among the models with
\eqn{1 \le r \le min(rank(X), jy*d)}. We then fix the local rank at
\eqn{\hat r} and repeat a similar procedure to determine \eqn{\hat rx}
and \eqn{\hat ry}, one at a time. Finally, with fixed \eqn{\hat rx} and \eqn{\hat ry},
we refine the estimation of r.
}
\examples{
library(NRRR)
set.seed(3)
# Simulation setting 2 in NRRR paper
simDat <- NRRR.sim(n = 100, ns = 100, nt = 100, r = 3, rx = 3, ry = 3,
                   jx = 8, jy = 8, p = 20, d = 20, s2n = 1, rho_X = 0.5,
                   rho_E = 0, Sigma = "CorrAR")
# using R function
fit_R <- with(simDat, NRRR.ic(Yest, Xest, Ag0 = NULL, Bg0 = NULL,
              jx = 8, jy = 8, p = 20, d = 20, n = 100,
              maxiter = 300, conv = 1e-4, method = c("RRR", "RRS")[1],
              lambda = 0, ic = c("BIC", "BICP", "AIC", "GCV")[1],
              dimred = c(TRUE, TRUE, TRUE), rankfix = NULL,
              xrankfix = NULL, yrankfix = NULL,
              lang=c('R','Rcpp')[1]))
# using Rcpp function
fit_Rcpp <- with(simDat, NRRR.ic(Yest, Xest, Ag0 = NULL, Bg0 = NULL,
                 jx = 8, jy = 8, p = 20, d = 20, n = 100,
                 maxiter = 300, conv = 1e-4, method = c("RRR", "RRS")[1],
                 lambda = 0, ic = c("BIC", "BICP", "AIC", "GCV")[1],
                 dimred = c(TRUE, TRUE, TRUE), rankfix = NULL,
                 xrankfix = NULL, yrankfix = NULL,
                 lang=c('R','Rcpp')[2]))
}
\references{
Liu, X., Ma, S., & Chen, K. (2020).
Multivariate Functional Regression via Nested Reduced-Rank Regularization.
arXiv: Methodology.
}
