#' Generate prediction of response curves
#'
#' This function generates prediction of response curves.
#'
#' @param sseq A sequence of the observed time points of x(s).
#' @param tseq A sequence of the prediction time points of y(t).
#' @param X X(s), an array of observed predictor curves.
#' @param phi The set of basis functions to expand x(s).
#' @param Ag i.e. U, a d*ry matrix of rank ry
#' @param Ag Global low-dimensional structure U, a d by ry matrix of rank ry.
#' @param Bg Global low-dimensional structure V, a p by rx matrix of rank rx.
#' @param Al Local low-dimensional structure A, a jy*ry by r matrix of rank r.
#' @param Bl Local low-dimensional structure B, a jx*rx by r matrix of rank r.
#' @return The returned items
#'   \item{Cstar}{\eqn{(U \otimes I_jy)A* B*^T(V \otimes I_jx)^T}, the coefficient matrix
#'           in equation (7).}
#'   \item{Jpsihalf}{psi is the set of basis functions to expand y(t); Jpsi is
#'           the covariance matrix of psi, and Jpsihalf is (Jpsi)^(1/2).}
#'   \item{Ypred}{the predicted response curves Ypred(t), an array of
#'           dimension (n, d, nt).}
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#' @importFrom splines bs
#' @export
#' @examples
#' library(NRRR)
#' set.seed(3)
#' # Simulation setting 2 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=100,nt=100,r=3,rx=3,ry=3,
#'                    jx=8,jy=8,p=20,d=20,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_nrrr <- with(simDat, NestRRR.select(Yest,Xest,Ag0=NULL,Bg0=NULL,
#'                               jx=8,jy=8,p=20,d=20,n=100,
#'                               maxiter=300,conv=1e-4,quietly=FALSE,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),rankfix=NULL))
#' Ypred_nrrr <- NestRRR.prediction(simDat$tseq,simDat$X,simDat$sseq,
#'                                  fit_nrrr$Ag,fit_nrrr$Bg,fit_nrrr$Al,fit_nrrr$Bl,
#'                                  simDat$phi)
NestRRR.prediction <- function(tseq,X,sseq,Ag,Bg,Al,Bl,phi){
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
  phi <- splines::bs(c(0,sseq),df = jx)[-1,]
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
