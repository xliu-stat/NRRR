% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RRR.r
\name{RRR}
\alias{RRR}
\title{Reduced-rank regression with a given rank}
\usage{
RRR(Y, X, nrank = 1, weight = FALSE, Gamma = diag(ncol(Y)),
    ypy.svd = TRUE, c.svd = FALSE)
}
\arguments{
\item{Y}{response matrix.}

\item{X}{design matrix.}

\item{nrank}{a user-specified rank. Default is 1.}

\item{weight}{a logical value. If TRUE, then the weighted criterion is employed.
Default is FALSE.}

\item{Gamma}{a weight matrix. Default is an identity matrix.}

\item{ypy.svd}{a logical value. If TRUE (default), then the singular value decomposition
method is used. If FALSE, then the eigenvalue decomposition method is used.}

\item{c.svd}{a logical value. If TRUE, then output the singular value
decomposition of the estimated coefficient
matrix. Default is FALSE.}
}
\value{
The function returns a list:
  \item{C}{the reduced-rank estimator.}
  \item{C_ls}{the least square estimator.}
  \item{rank}{the rank value.}
}
\description{
This function provides the reduced-rank regression estimator with a given rank.
}
\references{
Velu, R., & Reinsel, G. C. (2013). Multivariate reduced-rank regression:
theory and applications (Vol. 136). Springer Science & Business Media.
}
