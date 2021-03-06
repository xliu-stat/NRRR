% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NRRR.func.r
\name{NRRR.func}
\alias{NRRR.func}
\title{Multivatiate functional regression via nested reduced-rank regularization}
\usage{
NRRR.func(Y, X, jx, jy, degree = 3, S.interval = NULL, T.interval = NULL,
          tuning = c('CV', 'BIC', 'BICP', 'AIC', 'GCV')[1],
          nfold = 10, norder = NULL, method = c('RRR', 'RRS')[1],
          lambda = 0, maxiter = 300, conv = 1e-4,
          dimred = c(TRUE,TRUE,TRUE), rankfix = NULL, xrankfix = NULL,
          yrankfix = NULL, lang = c('R', 'Rcpp')[1])
}
\arguments{
\item{Y}{a data frame of the functional responses with d columns
for the values of the response components, a column indicating subject ID named as 'ID'
and a column of observation times named as 'TIME'.}

\item{X}{a data frame of the functional predictors with p columns
for the values of the predictor components, a column indicating subject ID named as 'ID'
and a column of observation times named as 'TIME'. Order of Subject ID should be the SAME as that of Y.}

\item{jx}{number of basis functions to expand the predictor trajectory.}

\item{jy}{number of basis functions to expand the response trajectory.}

\item{degree}{degree of piecewise polynomial. Default is 3. See \code{\link{bs}}.}

\item{S.interval}{range of observation times of X, e.g., \code{S.interval = c(0, 1)}. Default is NULL.}

\item{T.interval}{range of observation times of Y, e.g., \code{T.interval = c(0, 1)}. Default is NULL.}

\item{tuning}{methods to select ranks. If \code{tuning = 'CV'}, then cross validation
is used to select ranks. Otherwise, the selected information criterion
is used to select ranks.}

\item{nfold}{the number of folds used in cross validation. Default is 10.}

\item{norder}{a vector of length n that assigns samples to multiple folds for cross validation.}

\item{method}{'RRR' (default): no additional ridge penalty; 'RRS': add an
additional ridge penalty.}

\item{lambda}{the tuning parameter to control the amount of ridge
penalization. It is only used when \code{method = 'RRS'}.
Default is 0.}

\item{maxiter}{the maximum iteration number of the
blockwise coordinate descent algorithm. Default is 300.}

\item{conv}{the tolerance level used to control the convergence of the
blockwise coordinate descent algorithm. Default is 1e-4.}

\item{dimred}{a vector of logical values to decide whether to use the specified tuning method
do rank selection on certain dimensions. TRUE means the rank is selected
by the specified tuning method.
If \code{dimred[1] = FALSE}, r is provided by \code{rankfix}
or \eqn{min(jy*d, rank(X))};
If \code{dimred[2] = FALSE}, rx equals to \code{xrankfix} or p;
If \code{dimred[3] = FALSE},
ry equals to \code{yrankfix} or d.
Default is \code{c(TRUE, TRUE, TRUE)}.}

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
  \item{df}{the estimated degrees of freedom of the selected model.}
  \item{sse}{the sum of squared errors of the selected model.}
  \item{ic}{a vector containing values of BIC, BICP, AIC, GCV of the selected model.}
  \item{rx_path}{a matrix displays the path of selecting rx with cross validation. Only available when 'CV' is selected.}
  \item{ry_path}{a matrix displays the path of selecting ry with cross validation. Only available when 'CV' is selected.}
  \item{r_path}{a matrix displays the path of selecting r with cross validation. Only available when 'CV' is selected.}
  \item{rank}{the estimated r.}
  \item{rx}{the estimated rx.}
  \item{ry}{the estimated ry.}
  \item{sseq}{sequence of the time points of observing the predictor trajectory.}
  \item{phi}{the basis functions to expand the predictor trajectory.}
  \item{tseq}{sequence of the time points of observing the response trajectory.}
  \item{psi}{the basis functions to expand the response trajectory.}
}
\description{
This function takes functional observations as input and fits a nested
reduced-rank regression with a user-specified rank selection method.
The B-spline basis is used to conduct basis expansion.
}
\details{
This function applies a basis expansion and truncation method to transform
the functional problem into a classical finite-dimensional regression problem,
and then fits a nested reduced-rank regression.
The functional observations are commonly collected in a discrete mannar, and before using this
function, the functional response observations and functional predictor observations
should be saved as a data frame. Data standardization procedures, e.g., centering or
scaling, should be conducted before using this function. B-spline basis is
used to conduct basis expansion. If the functional data are already processed into an integrated form,
functions like \code{\link{NRRR.est}}, \code{\link{NRRR.ic}} or \code{\link{NRRR.cv}} can
be used to fit a NRRR model.
}
\examples{
n <- 100
ns <- 80
nt <- 80
p <- 10
d <- 10
jx <- 8
jy <- 8
library(NRRR)
set.seed(3)
# generate functional data
simDat <- NRRR.sim(n = 100, ns = 80, nt = 80, r = 3, rx = 3, ry = 3,
                   jx = 8, jy = 8, p = 10, d = 10, s2n = 1, rho_X = 0.5,
                   rho_E = 0, Sigma = "CorrAR")
# dimension: c(n, d, nt)
dim(simDat$Y)
# dimension: c(n, p, ns)
dim(simDat$X)
# process the functional data into a data frame with columns:
# predictor/response components, subject ID, and observation time

Y <- t(simDat$Y[1,,])
for (i in 2:n){
  Y <- rbind(Y, t(simDat$Y[i,,]))
}
Y <- data.frame(Y)
Y$TIME <- rep(simDat$tseq, n)
Y$ID <- rep(1:n, each = nt) # Y: nt*n rows, d+2 columns
head(Y)

X <- t(simDat$X[1,,])
for (i in 2:n){
  X <- rbind(X, t(simDat$X[i,,]))
}
X <- data.frame(X)
X$TIME <- rep(simDat$sseq, n)
X$ID <- rep(1:n, each = ns) # X: ns*n rows, p+2 columns
head(X)

fit_func <- NRRR.func(Y, X, jx, jy, degree = 3, S.interval = NULL, T.interval = NULL,
                      tuning = c('CV', 'BIC', 'BICP', 'AIC', 'GCV')[2],
                      nfold = 10, norder = NULL, method = c('RRR', 'RRS')[1],
                      lambda = 0, maxiter = 300, conv = 1e-4,
                      dimred = c(TRUE,TRUE,TRUE), rankfix = NULL, xrankfix = NULL,
                      yrankfix = NULL, lang = c('R', 'Rcpp')[1])
}
\references{
Liu, X., Ma, S., & Chen, K. (2020).
Multivariate Functional Regression via Nested Reduced-Rank Regularization.
arXiv: Methodology.
}
