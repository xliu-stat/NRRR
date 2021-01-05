#' @title
#' Predict the response trajectory
#'
#' @description
#' This function generates prediction of the functional response based on the fitting
#' results of the NRRR model and the given predictor trajectory. The functional
#' response is predicted at a given sequence of time points \code{tseq}.
#'
#' @usage
#' NRRR.pred(tseq, X, sseq, Ag, Bg, Al, Bl, phi)
#'
#' @param tseq a sequence of time points at which the prediction of the response
#'             trajectory is generated.
#' @param sseq a sequence of time points at which the predictor trajectory is observed.
#' @param X an array of dimension \code{(n, p, length(sseq))} where n is
#'          the sample size and p is the number of components in the
#'          multivariate predictor. It is the predictor trajectory observed at
#'          discrete time points \code{sseq}.
#' @param phi a matrix of dimension \code{length(sseq)}-by-jx. It is the set of
#'            basis functions to expand the predictor trajectory.
#' @param Ag the estimated matrix U in a NRRR model.
#' @param Bg the estimated matrix V in a NRRR model.
#' @param Al the estimated matrix A in a NRRR model.
#' @param Bl the estimated matrix B in a NRRR model.
#'
#'
#' @return This function returns a list:
#'   \item{Cstar}{the estimated coefficient matrix in equation (7) of the NRRR paper,
#'                i.e., \eqn{(U \otimes I_jy)A* B*^T(V \otimes I_jx)^T}.}
#'   \item{Ypred}{the predicted response trajectory which is an array of
#'           dimension \code{(n, d, length(tseq))}}
#'
#'
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#'
#' @importFrom splines bs
#' @export
#' @examples
#' library(NRRR)
#' set.seed(3)
#' # Simulation setting 2 in NRRR paper
#' simDat <- NRRR.sim(n = 100, ns = 100, nt = 100, r = 3, rx = 3, ry = 3,
#'                    jx = 8, jy = 8, p = 20, d = 20, s2n = 1, rho_X = 0.5,
#'                    rho_E = 0, Sigma = "CorrAR")
#' fit_nrrr <- with(simDat, NRRR.ic(Yest, Xest, Ag0 = NULL, Bg0 = NULL,
#'                  jx = 8, jy = 8, p = 20, d = 20, n = 100,
#'                  maxiter = 300, conv = 1e-4, method = c("RRR", "RRS")[1],
#'                  lambda = 0, ic = c("BIC", "BICP", "AIC", "GCV")[1],
#'                  dimred = c(TRUE, TRUE, TRUE), rankfix = NULL))
#' Ypred_nrrr <- NRRR.pred(simDat$tseq, simDat$X, simDat$sseq,
#'                         fit_nrrr$Ag, fit_nrrr$Bg, fit_nrrr$Al,
#'                         fit_nrrr$Bl, simDat$phi)
NRRR.pred <- function(tseq,X,sseq,Ag,Bg,Al,Bl,phi){
  n <- nrow(X)

  p <- nrow(Bg)
  d <- nrow(Ag)
  rx <- ncol(Bg)
  ry <- ncol(Ag)

  jx <- nrow(Bl)/rx
  jy <- nrow(Al)/ry

  nt <- length(tseq)
  ns <- length(sseq)

  # Note that you can predict at different time
  psi <- splines::bs(c(0,tseq),df = jy)[-1,]
  Jpsi <- matrix(nrow=jy,ncol=jy,0)
  tdiff <- (tseq - c(0,tseq[-nt]))
  for(t in 1:nt) Jpsi <- Jpsi + psi[t,]%*%t(psi[t,])*tdiff[t]
  eJpsi <- eigen(Jpsi)
  Jpsihalf <- eJpsi$vectors%*%diag(sqrt(eJpsi$values))%*%t(eJpsi$vectors)
  Jpsihalfinv <- eJpsi$vectors%*%diag(1/sqrt(eJpsi$values))%*%t(eJpsi$vectors)

  # Compute Cstar from estimated C
  alindex <- rep(1:ry,jy)
  Alstar <- Al[order(alindex),]
  blindex <- rep(1:rx,jx)
  Blstar <- Bl[order(blindex),]
  # Adjust
  Alstar <- kronecker(diag(ry),Jpsihalfinv)%*%Alstar
  Cstar <- kronecker(Ag,diag(jy))%*%Alstar%*%t(Blstar)%*%kronecker(t(Bg),diag(jx))

  # Integrated X
  # phi <- splines::bs(c(0,sseq),df = jx)[-1,]
  Xint <- matrix(nrow=p*jx,ncol=n,0)
  sdiff <- (sseq - c(0,sseq[-ns]))
  for(s in 1:ns){
    Xint <- Xint + kronecker(diag(nrow=p,ncol=p),phi[s,])%*%t(X[,,s])*sdiff[s]
  }


  # Compute predicted values
  Ypred <- array(dim=c(n,d,nt),NA)
  for(t in 1:nt){
    Psit <- kronecker(diag(nrow=d,ncol=d),t(psi[t,]))
    Ypred[,,t] <- t(Psit%*%Cstar%*%Xint)
  }

  list(Jpsihalf=Jpsihalf,Cstar=Cstar,Ypred=Ypred)

}
