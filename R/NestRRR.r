#' Obtain NRRR estimator with given rank
#'
#' This function uses the proposed blockwise coordinate descent algorithm
#'  to get the NRRR estimator with given rank (r, rx, ry).
#'
#' @param Y Response matrix with n rows and jy*d columns.
#' @param X Design matrix with n rows and jx*p columns.
#' @param Ag0 Initial estimator of U, if NULL then generate it by
#'           function NestRRRini.
#' @param Bg0 Initial estimator of V, if NULL then generate it by
#'           function NestRRRini.
#' @param rini,r Rank of the local reduced-rank structure. \emph{rini} is used
#'               in \code{NestRRRinit()} to get the initial estimator of U and V.
#' @param rx Number of latent predictors.
#' @param ry Number of latent responses.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param n Sample size.
#' @param maxiter Maximum iteration number, default 100.
#' @param conv Tolerance level to control convergence, default 1E-4.
#' @param quietly FALSE show the final fitting information (SSE, BIC, df), if TRUE (default) then not show.
#' @param method RRR (default) no ridge penalty; if RRS then use ridge penalty.
#' @param lambda Tuning parameter for the ridge penalty, only used when method='RRS', default=0.
#' @return The returned items
#'   \item{Ag}{the global low-dimensional structure U, a d by ry matrix of rank ry.}
#'   \item{Bg}{the global low-dimensional structure V, a p by rx matrix of rank rx.}
#'   \item{Al}{the local low-dimensional structure A.}
#'   \item{Bl}{the local low-dimensional structure B.}
#'   \item{C}{the NRRR estimator of the coefficient matrix C.}
#'   \item{df}{a scalar, the estimated degrees of freedom of the NRRR model.}
#'   \item{sse}{a scalar, the objective function value.}
#'   \item{ic}{a vector contains values of BIC,BICP,AIC,GCV.}
#'   \item{obj}{ a vector contains all objective function (sse) values along iteration.}
#'   \item{iter}{a scalar, number of iterations to converge.}
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#' @examples
#' library(NRRR)
#' simDat <- nrrr.sim(
#'   n = 100, ns = 200, nt = 200, r = 5, rx = 3, ry = 3,
#'   jx = 15, jy = 15, p = 10, d = 6, s2n = 1, rho_X = 0.5,
#'   rho_E = 0, Sigma = CorrAR
#' )
#' fit_init <- with(simDat, NestRRR(
#'   Y = Yest, X = Xest, Ag0 = NULL, Bg0 = NULL,
#'   rini = 5, r = 5,
#'   rx = 3, ry = 3, jx = 15, jy = 15, p = 10, d = 6, n = 100
#' ))
#' fit_init$Ag
#' @importFrom MASS ginv
#' @export
NestRRR <- function(Y, X, Ag0 = NULL, Bg0 = NULL, rini, r, rx, ry, jx, jy, p, d, n,
                    maxiter = 100, conv = 1e-4, quietly = TRUE,
                    method = c("RRR", "RRS")[1], lambda = 0) {
  if (method == "RRS" & lambda != 0) { # RRS is reduced rank ridge regression.
    Y <- rbind(Y, matrix(nrow = p * jx, ncol = d * jy, 0.0))
    X <- rbind(X, sqrt(lambda) * diag(nrow = p * jx, ncol = p * jx))
    n <- n + p * jx
  }

  # require(rrpack)
  # compute initial values of U and V
  if (is.null(Ag0) | is.null(Bg0)) {
    ini <- NestRRRini(Y, X, rini, rx, ry, jx, jy, p, d, n)
    Ag0 <- ini$Ag
    Bg0 <- ini$Bg
  }

  if (d == ry) {
    Yl0 <- Y
  } else {
    Yl0 <- Y %*% kronecker(diag(jy), Ag0)
  }
  if (p == rx) {
    Xl0 <- X
  } else {
    Xl0 <- X %*% kronecker(diag(jx), Bg0)
  }

  # Given Ag (U) and Bg (V), compute Al and Bl
  fitRR <- RRR(Yl0, Xl0, nrank = r)
  Bl0 <- fitRR$C_ls %*% fitRR$A
  Al0 <- fitRR$A


  obj <- vector() # collect all objective function (sse) values
  obj[1] <- Obj(Y, X, Ag0, Bg0, Al0, Bl0, jx, jy)$obj
  objnow <- obj[1] + 10
  iter <- 1
  while (iter < maxiter & abs(obj[iter] - objnow) > conv) { # use sse to control convergence.

    ### For updating Ag (U)###
    if (d == ry) {
      Ag1 <- diag(nrow = ry, ncol = ry) # an identity matrix
      Yl0 <- Y
    } else {
      Xg <- Xl0 %*% Bl0 %*% t(Al0)
      Yg0 <- matrix(nrow = d, ncol = ry, 0)
      for (j in 1:jy) {
        a1 <- d * (j - 1) + 1
        b1 <- d * j

        a2 <- ry * (j - 1) + 1
        b2 <- ry * j

        Yg0 <- Yg0 + t(Y[, a1:b1]) %*% (Xg[, a2:b2])
      }
      svdYg0 <- svd(Yg0)
      Ag1 <- svdYg0$u %*% t(svdYg0$v)
      Yl0 <- Y %*% kronecker(diag(jy), Ag1)
    }


    ### For updating Bg (V)###
    if (p == rx) {
      Bg1 <- diag(nrow = rx, ncol = rx)
      Xl0 <- X
    } else {
      yB <- as.vector(Yl0 %*% Al0)
      XB <- matrix(nrow = n * r, ncol = rx * p, 0)
      for (j in 1:jx) {
        a1 <- rx * (j - 1) + 1
        b1 <- rx * j

        a2 <- p * (j - 1) + 1
        b2 <- p * j

        XB <- XB + kronecker(matrix(t(Bl0)[, a1:b1], nrow = r, ncol = rx), X[, a2:b2])
      }
      Bg1 <- matrix(nrow = p, ncol = rx, byrow = FALSE, MASS::ginv(t(XB) %*% XB +
                   0.0000001 * diag(1, rx * p, rx * p)) %*% t(XB) %*% yB)

      Bg1 <- qr.Q(qr(Bg1))
      Xl0 <- X %*% kronecker(diag(jx), Bg1)
    }

    ### For updating Al (A) and Bl (B)###
    fitRR <- RRR(Yl0, Xl0, nrank = r)
    Bl1 <- fitRR$C_ls %*% fitRR$A
    Al1 <- fitRR$A

    iter <- iter + 1
    Objiter <- Obj(Y, X, Ag1, Bg1, Al1, Bl1, jx, jy)
    obj[iter] <- Objiter$obj
    objnow <- obj[iter - 1]
    C <- Objiter$C


    Al0 <- Al1
    Bl0 <- Bl1
    Ag0 <- Ag1
    Bg0 <- Bg1
  }

  # compute rank of X
  xr <- sum(svd(X)$d > 1e-2)
  # compute degrees of freedom of NRRR
  df <- ifelse(xr / jx > rx, rx * (xr / jx - rx), 0) + ry * (d - ry) + (jy * ry + jx * rx - r) * r
  # obtain sse
  sse <- ifelse(Objiter$sse < 0.1, 0, Objiter$sse)
  # compute information criterion
  BIC <- log(sse) + log(d * jy * n) / d / jy / n * df
  BICP <- log(sse) + 2 * log(d * jy * p * jx) / d / jy / n * df
  AIC <- log(sse) + 2 / d / jy / n * df
  GCV <- sse / d / jy / n / (1 - df / d / jy / n)^2

  if (method == "RRS") sse <- sse - lambda * sum(Objiter$C^2)

  if (!quietly) {
    cat("SSE = ", sse, "   BIC = ", BIC, "   DF = ", df, "\n", sep = "")
  }
  return(list(
    Ag = Ag1, Bg = Bg1, Al = as.matrix(Al1), Bl = as.matrix(Bl1), C = C, df = df,
    sse = sse, ic = c(BIC, BICP, AIC, GCV), obj = obj, iter = iter
  ))
}
