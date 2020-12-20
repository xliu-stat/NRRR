#' @title
#' Nested reduced-rank regression with a given rank
#'
#' @description
#' This function uses a blockwise coordinate descent algorithm
#' to get the nested reduced-rank regression estimator with a
#' given \code{(r, rx, ry)}.
#'
#' @usage
#' NestRRR(Y, X, Ag0 = NULL, Bg0 = NULL, rini, r, rx, ry, jx, jy, p, d, n,
#'         maxiter = 100, conv = 1e-4, quietly = TRUE,
#'         method = c("RRR", "RRS")[1], lambda = 0)
#'
#'
#' @param Y the response matrix of dimension n-by-jy*d.
#' @param X the design matrix of dimension n-by-jx*p.
#' @param Ag0 an initial estimator of matrix U. If NULL then generate it
#'            by \code{\link{NestRRRini}}. Default is NULL.
#' @param Bg0 an initial estimator of matrix V, if NULL then generate it
#'            by \code{\link{NestRRRini}}. Default is NULL.
#' @param rini,r rank of the local reduced-rank structure. \code{rini} is used
#'               in \code{\link{NestRRRini}} to get the initial
#'               estimator of U and V.
#' @param rx the number of latent predictors.
#' @param ry the number of latent responses.
#' @param jx the number of basis functions to expand functional predictor.
#' @param jy the number of basis functions to expand functional response.
#' @param p the number of predictors.
#' @param d the number of responses.
#' @param n the sample size.
#' @param maxiter the maximum iteration number of the
#'                blockwise coordinate descent algorithm. Default is 100.
#' @param conv the tolerance level used to control the convergence of the
#'             blockwise coordinate descent algorithm. Default is 1e-4.
#' @param quietly a logical value. FALSE: show the
#'                final fitting information (sse and df);
#'                TRUE (default): do not show the results.
#' @param method 'RRR' (default): no additional ridge penalty; 'RRS': add an
#'               additional ridge penalty.
#' @param lambda the tuning parameter to control the amount of ridge
#'               penalization. It is only used when \code{method = 'RRS'}.
#'               Default is 0.
#'
#' @return The function returns a list:
#'   \item{Ag}{the estimated U.}
#'   \item{Bg}{the estimated V.}
#'   \item{Al}{the estimated A.}
#'   \item{Bl}{the estimated B.}
#'   \item{C}{the estimated coefficient matrix C.}
#'   \item{df}{degrees of freedom of the model.}
#'   \item{sse}{the sum of squared errors.}
#'   \item{ic}{a vector containing values of BIC, BICP, AIC, GCV.}
#   \item{obj}{a vector contains all objective function (sse) values along iteration.}
#'   \item{iter}{the number of iterations to converge.}
#'
#' @details
#' The \emph{nested reduced-rank regression (NRRR)} is first motivated to solve
#' a multivariate functional linear regression problem where both the response and
#' predictor are multivariate and functional (i.e., \eqn{Y(t)=(y_1(t),...,y_d(t))^T}
#' and \eqn{X(s)=(x_1(s),...,y_p(s))^T}). To control the complexity of the
#' problem, NRRR proposes a nested reduced-rank structure on the regression
#' surface \eqn{C(s,t)}. Specifically, a global dimension reduction makes use of the
#' correlation within the components of multivariate response and multivariate
#' predictor. Matrices U (d-by-ry) and V (p-by-rx) provides weights to form latent
#' functional responses and latent functional predictors. Dimension reduction is achieved
#' once \eqn{ry \le d} or \eqn{rx \le p}. Then, a local dimension reduction is
#' conducted by restricting the latent regression surface \eqn{C^*(s,t)} to be of low-rank.
#' After basis expansion and truncation, also by applying proper rearrangement to
#' columns and rows of the resulting data matrices and coefficient matrices, we
#' have the nested reduced-rank problem:
#' \deqn{ \min_{C} || Y - XC ||_F^2, s.t., C = (I_{jx} \otimes V) BA^T (I_{jy} \otimes U)^T,}
#' where \eqn{BA^T} is a full-rank decomposition structure to control the local
#' rank and \eqn{jx, jy} are the number of basis functions. Beyond the functional
#' setup, this structure can also be applied in multiple scenarios, including
#' multivariate time series autoregression analysis and tensor-on-tensor regression.
#' This problem is non-convex and has no explicit solution, thus we use a
#' blockwise coordinate descent algorithm to find a local solution.
#'
#'
#'
#'
#'
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
    cat("SSE = ", sse, "DF = ", df, "\n", sep = "")
  }
  return(list(
    Ag = Ag1, Bg = Bg1, Al = as.matrix(Al1), Bl = as.matrix(Bl1), C = C, df = df,
    sse = sse, ic = c(BIC, BICP, AIC, GCV), obj = obj, iter = iter
  ))
}
