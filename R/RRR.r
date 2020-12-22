#' @title
#' Reduced-rank regression with a given rank
#'
#' @description
#' This function provides the reduced-rank regression estimator with a given rank.
#'
#' @usage
#' RRR(Y, X, nrank = 1, weight = FALSE, Gamma = diag(ncol(Y)),
#'     ypy.svd = TRUE, c.svd = FALSE)
#'
#' @param Y response matrix.
#' @param X design matrix.
#' @param nrank a user-specified rank. Default is 1.
#' @param weight a logical value. If TRUE, then weighted criterion is performed.
#'               Default is FALSE.
#' @param Gamma a weight matrix. Default is an identity matrix.
#' @param ypy.svd a logical value. If TRUE, svd function is used. Default is FALSE.
#' @param c.svd a logical value. If TRUE, output the svd of the coefficient
#'              matrix. Default is FALSE.
#'
#' @return The function returns a list:
#'   \item{C}{the reduced-rank estimator.}
#'   \item{C_ls}{the least square estimator.}
#'   \item{rank}{the rank value.}
#'
#' @references
#' Velu, R., & Reinsel, G. C. (2013). Multivariate reduced-rank regression:
#' theory and applications (Vol. 136). Springer Science & Business Media.
#' @importFrom MASS ginv
#' @export
RRR <- function(Y, X, nrank = 1, weight = FALSE, Gamma = diag(ncol(Y)),
                 ypy.svd = TRUE, c.svd = FALSE){
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  S_yx <- t(Y) %*% X
  S_xx <- t(X) %*% X
  S_xx_inv <- tryCatch(MASS::ginv(S_xx), error = function(e) solve(S_xx +
                                                               0.1 * diag(p)))
  if (sum(is.na(S_xx_inv)) > 0) {
    S_xx_inv <- solve(S_xx + 0.1 * diag(p))
  }
  C_ls <- S_xx_inv %*% t(S_yx)
  if (weight == TRUE) {
    eigenGm <- eigen(Gamma)
    sqrtGm <- eigenGm$vectors %*% diag(sqrt(eigenGm$values)) %*%
      t(eigenGm$vectors)
    sqrtinvGm <- eigenGm$vectors %*% diag(sqrt(eigenGm$values^(-1))) %*%
      t(eigenGm$vectors)
    if (ypy.svd) {
      XC <- X %*% C_ls %*% sqrtGm
      svdXC <- tryCatch(svd(XC, nu = nrank, nv = nrank),
                        error = function(e) 2)
      if (tryCatch(svdXC == 2, error = function(e) 3) ==
          3) {
        A <- svdXC$v[, 1:nrank]
        Ad <- (svdXC$d[1:nrank])^2
      }
      else {
        ypy.svd <- FALSE
      }
    }
    if (!ypy.svd) {
      SS <- sqrtGm %*% S_yx %*% C_ls %*% sqrtGm
      SS <- (SS + t(SS))/2
      eigenSS <- eigen(SS, symmetric = TRUE)
      A <- as.matrix(eigenSS$vectors[, 1:nrank])
      Ad <- eigenSS$values[1:nrank]
    }
    AA <- A %*% t(A)
    C_rr <- C_ls %*% sqrtGm %*% AA %*% sqrtinvGm
  }
  else {
    if (ypy.svd) {
      XC <- X %*% C_ls
      svdXC <- tryCatch(svd(XC, nu = nrank, nv = nrank),
                        error = function(e) 2)
      if (tryCatch(svdXC == 2, error = function(e) 3) ==
          3) {
        A <- svdXC$v[, 1:nrank]
        Ad <- (svdXC$d[1:nrank])^2
      }
      else {
        ypy.svd <- FALSE
      }
    }
    if (!ypy.svd) {
      SS <- S_yx %*% C_ls
      SS <- (SS + t(SS))/2
      eigenSS <- eigen(SS, symmetric = TRUE)
      A <- as.matrix(eigenSS$vectors[, 1:nrank])
      Ad <- eigenSS$values[1:nrank]
    }
    AA <- A %*% t(A)
    C_rr <- C_ls %*% AA
  }
  if (c.svd) {
    svd_C <- svd(C_rr, nv = nrank, nu = nrank)
    U <- as.matrix(svd_C$u[, 1:nrank])
    V <- as.matrix(svd_C$v[, 1:nrank])
    D <- diag(svd_C$d[1:nrank], nrow = nrank)
    list(A = A, Ad = Ad, C_ls = C_ls, C_rr = C_rr, U = U,
         V = V, D = D, C = C_rr, rank = nrank)
  }
  else {
    list(A = A, Ad = Ad, C_ls = C_ls, C_rr = C_rr, C = C_rr,
         rank = nrank)
  }
}
