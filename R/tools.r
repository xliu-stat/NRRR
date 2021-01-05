# Covariance structure: Compound symmetry
#
# This function generates a covariance matrix with compound symmetry structure.
#
# @param p Dimension.
# @param rho Correlation strength.
# @return A covariance matrix (p by p).
# @export
# @examples
# library(NRRR)
# CorrCS(10, 0.5)
CorrCS <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,rho)
  diag(Sigma) <- 1
  Sigma
}

# Covariance structure: Autoregressive
#
# This function generates a covariance matrix with autoregressive structure.
#
# @param p Dimension.
# @param rho Correlation strength.
# @return A covariance matrix (p by p).
# @export
# @examples
# library(NRRR)
# CorrAR(10, 0.5)
CorrAR <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,NA)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  Sigma
}

# Compute || Y - XC ||_F^2
#
# This function computes the sum of squared errors
# with given (A,B,U,V) and (X,Y).
#
# @param Y Response matrix with n rows and jy*d columns.
# @param X Design matrix with n rows and jx*p columns.
# @param Ag Matrix U.
# @param Bg Matrix V.
# @param Al Matrix A.
# @param Bl Matrix B.
# @param jx Number of basis functions to expand x(s).
# @param jy Number of basis functions to expand y(t).
# @export
# @return The returned items
#   \item{C}{the NRRR estimator of the coefficient matrix C.}
#   \item{XC}{a matrix, prediction of Y.}
#   \item{E}{a matrix, estimated error matrix.}
#   \item{sse,obj}{ a scaler, the sum of squared errors.}
Obj <- function(Y,X,Ag,Bg,Al,Bl,jx,jy){

  if(nrow(Bg)!=ncol(Bg)){
    Bl <- kronecker(diag(jx),Bg)%*%Bl
  }
  if(nrow(Ag)!=ncol(Ag)){
    Al <- kronecker(diag(jy),Ag)%*%Al
  }
  C <- Bl%*%t(Al)
  XC <- X%*%C
  E <- Y - XC
  sse <- sum(E^2)

  return(list(C=C,XC=XC,E=E,sse=sse,obj=sse))
}


cv.rrr <- function (Y,
                    X,
                    nfold = 10,
                    maxrank = min(dim(Y), dim(X)),
                    norder = NULL,
                    coefSVD = FALSE)
{
  ## record function call
  Call <- match.call()

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  ndel <- round(n / nfold)
  if (is.null(norder))
    norder <- sample(seq_len(n), n)

  cr_path <- matrix(ncol = nfold, nrow = maxrank + 1, NA)
  for (f in seq_len(nfold)) {
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    Xf <- X[-iddel,]
    Xfdel <- X[iddel,]
    Yf <- Y[-iddel,]
    Yfdel <- Y[iddel,]
    ini <- rrr.fit(Yf, Xf, nrank = maxrank, coefSVD =  coefSVD)
    C_ls <- ini$coef.ls
    A <- ini$A
    tempFit <- Xfdel %*% C_ls
    tempC <- matrix(nrow = q, ncol = q, 0)
    for (i in seq_len(maxrank)) {
      tempC <- tempC + tcrossprod(A[, i])
      cr_path[i + 1, f] <-
        sum((Yfdel - tempFit %*% tempC) ^ 2)
    }
    cr_path[1, f] <- sum(Yfdel ^ 2)
  }
  index <- order(colSums(cr_path))
  crerr <- rowSums(cr_path[, index]) / length(index) * nfold
  minid <- which.min(crerr)
  rankest <- minid - 1
  ini <- rrr.fit(Y, X, nrank = maxrank)
  C_ls <- ini$coef.ls
  A <- ini$A

  out <- if (identical(rankest, 0)) {
    list(
      call = Call,
      cr.path = cr_path,
      cr.error = crerr,
      norder = norder,
      coef = matrix(0, nrow = p, ncol = q),
      rank = 0,
      coef.ls = C_ls
    )
  }
  else {
    list(
      call = Call,
      cr.path = cr_path,
      cr.error = crerr,
      norder = norder,
      coef.ls = C_ls,
      coef = C_ls %*% tcrossprod(A[, seq_len(rankest)]),
      rank = rankest
    )
  }
  #class(out) <- "cv.rrr"
  out
}



rrr.fit <- function(Y,
                    X,
                    nrank = 1,
                    weight = NULL,
                    coefSVD = FALSE)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  stopifnot(n == nrow(X))

  S_yx <- crossprod(Y, X)
  S_xx <- crossprod(X)

  ## FIXME: 0.1 is too arbitrary
  S_xx_inv <- tryCatch(
    MASS::ginv(S_xx),
    error = function(e)
      solve(S_xx + 0.1 * diag(p))
  )

  ## FIXME: if weighted, this needs to be weighted too
  C_ls <- tcrossprod(S_xx_inv, S_yx)

  if (!is.null(weight)) {
    stopifnot(nrow(weight) == q && ncol(weight) == q)
    eigenGm <- eigen(weight)
    ## FIXME: ensure eigen success?
    ## sqrtGm <- tcrossprod(eigenGm$vectors * sqrt(eigenGm$values),
    ##                      eigenGm$vectors)
    ## sqrtinvGm <- tcrossprod(eigenGm$vectors / sqrt(eigenGm$values),
    ##                         eigenGm$vectors)
    sqrtGm <- eigenGm$vectors %*% (sqrt(eigenGm$values) *
                                     t(eigenGm$vectors))
    sqrtinvGm <- eigenGm$vectors %*% (1 / sqrt(eigenGm$values) *
                                        t(eigenGm$vectors))

    XC <- X %*% C_ls %*% sqrtGm
    ## FIXME: SVD may not converge
    ## svdXC <- tryCatch(svd(XC,nu=nrank,nv=nrank),error=function(e)2)
    svdXC <- svd(XC, nrank, nrank)
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
    AA <- tcrossprod(A)
    C_rr <- C_ls %*% sqrtGm %*% AA %*% sqrtinvGm
  } else {
    ## unweighted
    XC <- X %*% C_ls
    svdXC <- svd(XC, nrank, nrank)
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
    AA <- tcrossprod(A)
    C_rr <- C_ls %*% AA
  }

  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    rank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
    coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  #class(ret) <- "rrr.fit"
  ret
}

rrs.fit <- function(Y,
                    X,
                    nrank = min(ncol(Y), ncol(X)),
                    lambda = 1,
                    coefSVD = FALSE)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  S_yx <- t(Y) %*% X
  ## This is a key difference
  S_xx <- t(X) %*% X + lambda * diag(p)

  ## S_xx_inv <- tryCatch(solve(S_xx+0.1*diag(p)),error=function(e)ginv(S_xx))
  ## S_xx_inv <- ginv(S_xx)
  ## Use the Woodbury matrix identity
  if (lambda != 0) {
    S_xx_inv <- 1 / lambda * diag(p) -
      lambda ^ (-2) * t(X) %*% MASS::ginv(diag(n) + lambda ^ (-1) * X %*%
                                            t(X)) %*% X
  } else{
    S_xx_inv <- MASS::ginv(S_xx)
    if (sum(is.na(S_xx_inv)) > 0) {
      S_xx_inv <- solve(S_xx + 0.1 * diag(p))
    }
  }

  C_ls <- S_xx_inv %*% t(S_yx)

  ypy.svd <- TRUE
  ##if(ypy.svd){
  ##This is another key difference
  XC <- rbind(X, sqrt(lambda) * diag(p)) %*% C_ls
  svdXC <- tryCatch(
    svd(XC, nu = nrank, nv = nrank),
    error = function(e)
      2)
  if (tryCatch(
    svdXC == 2,
    error = function(e)
      3) == 3) {
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
  } else{
    ypy.svd <- FALSE
  }
  #}
  if (!ypy.svd) {
    SS <- S_yx %*% C_ls
    SS <- (SS + t(SS)) / 2
    eigenSS <- eigen(SS, symmetric = TRUE)
    A <- as.matrix(eigenSS$vectors[, 1:nrank])
    Ad <- eigenSS$values[1:nrank]
  }

  AA <- A %*% t(A)
  C_rr <- C_ls %*% AA

  ##    if(c.svd){
  ##      svd_C <- svd(C_rr,nv=nrank,nu=nrank)
  ##      U <- as.matrix(svd_C$u[,1:nrank])
  ##      V <- as.matrix(svd_C$v[,1:nrank])
  ##      D <- diag(svd_C$d[1:nrank],nrow=nrank)
  ##
  ##      ####return ls estimator C_ls, reduced-rank estimator C_rr
  ##      ####return SVD of C_rr
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,U=U,V=V,D=D,C=C_rr,rank=nrank)
  ##    }else{
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,C=C_rr,rank=nrank)
  ##    }

  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    nrank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  #class(ret) <- "rrs.fit"
  ret
}





