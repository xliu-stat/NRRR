#' @title
#' Reduced-rank ridge regression with rank and tuning parameter selected by cross validation
#'
#' @description
#' This function performs reduced-rank ridge regression with the rank and the
#' tuning parameter selected by cross validation.
#'
#' @usage
#' RRRR.cv(Y, X, nfold = 10, rankmax = min(dim(Y), dim(X)), nlam = 100,
#'         lambda = seq(0, 100, length = nlam), norder = NULL,
#'         nest.tune = FALSE, fold.drop = 0)
#'
#' @param Y response matrix.
#' @param X design matrix.
#' @param nfold the number of folds used in cross validation. Default is 10.
#' @param rankmax the maximum rank allowed.
#' @param nlam the number of tuning parameter candidates. Default is 100.
#' @param lambda the tuning sequence of length \code{nlam}.
#' @param norder a vector of length n that assigns samples to multiple folds
#'               for cross validation. Default is NULL and then \code{norder} will be
#'               generated randomly.
#' @param nest.tune a logical value to specify whether to tune the rank and lambda
#'                  in a nested way. Default is FALSE.
#' @param fold.drop the number of folds to drop. Default is 0.
#'
#'
#' @return The function returns a list:
#'   \item{cr_path}{a matrix displays the path of model selection.}
#'   \item{C}{the estimated low-rank coefficient matrix.}
#'   \item{rank}{the selected rank.}
#'   \item{lam}{the selected tuning parameter for ridge penalty.}
#'
#' @references Mukherjee, A., & Zhu, J. (2011).
#' Reduced Rank Ridge Regression and Its Kernel Extensions.
#' Statistical analysis and data mining,
#' 4(6), 612–622.
#'
#' @importFrom rrpack rrs.fit
#' @export
RRRR.cv<-function(Y, X, nfold = 10, rankmax = min(dim(Y), dim(X)), nlam = 100,
                     lambda = seq(0, 100, length = nlam), norder = NULL, nest.tune = FALSE,
                     fold.drop = 0)
{
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)
  ndel <- round(n/nfold)
  if (is.null(norder)) {
    norder <- sample(1:n, n)
  }
  cr_path <- array(dim = c(nlam, rankmax + 1, nfold), NA)
  for (f in 1:nfold) {
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    Xf <- X[-iddel, ]
    Xfdel <- X[iddel, ]
    Yf <- Y[-iddel, ]
    Yfdel <- Y[iddel, ]
    for (ll in 1:nlam) {
      ini <- rrpack::rrs.fit(Yf, Xf, lambda = lambda[ll], nrank = rankmax)
      C_ls <- ini$coef.ls
      A <- ini$A
      tempFit <- Xfdel %*% C_ls
      tempC <- matrix(nrow = nrow(A), ncol = nrow(A), 0)
      for (i in 1:rankmax) {
        tempC <- tempC + A[, i] %*% t(A[, i])
        cr_path[ll, i + 1, f] <- sum((Yfdel - tempFit %*%
                                        tempC)^2)
      }
      cr_path[ll, 1, f] <- sum(Yfdel^2)
    }
  }
  index <- order(apply(cr_path, 3, sum))[1:ifelse(nfold > 5,
                                                  nfold - fold.drop, nfold)]
  crerr <- apply(cr_path[, , index], c(1, 2), sum)/length(index) *
    nfold
  minid <- which.min(crerr)
  minid2 <- ceiling(minid/nlam)
  rankest <- minid2 - 1
  minid1 <- minid - nlam * (minid2 - 1)
  if (nest.tune == TRUE) {
    lambda <- exp(seq(log(lambda[max(minid1 - 1, 1)] + 1e-04),
                      log(lambda[min(minid1 + 1, nlam)]), length = nlam))
    cr_path <- array(dim = c(nlam, rankmax + 1, nfold), NA)
    ndel <- round(n/nfold)
    for (f in 1:nfold) {
      if (f != nfold) {
        iddel <- norder[(1 + ndel * (f - 1)):(ndel *
                                                f)]
      }
      else {
        iddel <- norder[(1 + ndel * (f - 1)):n]
        ndel <- length(iddel)
      }
      nf <- n - ndel
      Xf <- X[-iddel, ]
      Xfdel <- X[iddel, ]
      Yf <- Y[-iddel, ]
      Yfdel <- Y[iddel, ]
      for (ll in 1:nlam) {
        ini <- rrpack::rrs.fit(Yf, Xf, lambda = lambda[ll], nrank = rankmax)
        C_ls <- ini$coef.ls
        A <- ini$A
        tempFit <- Xfdel %*% C_ls
        tempC <- matrix(nrow = nrow(A), ncol = nrow(A),
                        0)
        for (i in 1:rankmax) {
          tempC <- tempC + A[, i] %*% t(A[, i])
          cr_path[ll, i + 1, f] <- sum((Yfdel - tempFit %*%
                                          tempC)^2)
        }
        cr_path[ll, 1, f] <- sum(Yfdel^2)
      }
    }
    index <- order(apply(cr_path, 3, sum))[1:ifelse(nfold >
                                                      5, nfold - fold.drop, nfold)]
    crerr <- apply(cr_path[, , index], c(1, 2), sum)/length(index) *
      nfold
    minid <- which.min(crerr)
    minid2 <- ceiling(minid/nlam)
    rankest <- minid2 - 1
    minid1 <- minid - nlam * (minid2 - 1)
  }
  if (rankest == 0) {
    list(cr_path = cr_path, CRE = crerr, norder = norder,
         C = matrix(nrow = p, ncol = q, 0), rank = 0)
  }
  else {
    fit <- rrpack::rrs.fit(Y, X, lambda = lambda[minid1], nrank = rankmax)
    C_lslam <- fit$coef.ls
    A <- fit$A
    list(cr_path = cr_path, CRE = crerr, ID = c(minid1, minid2),
         norder = norder, C = C_lslam %*% A[, 1:rankest] %*%
           t(A[, 1:rankest]), rank = rankest, lam = lambda[minid1])
  }
}
