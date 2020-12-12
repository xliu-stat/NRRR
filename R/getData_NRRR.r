#' Generate simulation data
#'
#' This function generates simulation data for the nested reduced-rank
#' functional regression model.
#'
#' @param n Sample size.
#' @param ns Number of discrete observations of x(s).
#' @param nt Number of discrete observations of y(t).
#' @param r Rank.
#' @param rx Number of latent predictors.
#' @param ry Number of latent responses.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param s2n Signal-to-noise ratio.
#' @param rho_X The correlation strength among columns of X,
#'              a scalar between 0 and 1.
#' @param rho_E The correlation strength among columns of E,
#'              a scalar between 0 and 1.
#' @param Sigma The correlation structure, select from
#'           CorrAR and CorrCS.
#' @return The generated data containing
#'   \item{Ag}{the global low-dimensional structure U, a d by ry matrix of rank ry.}
#'   \item{Bg}{the global low-dimensional structure V, a p by rx matrix of rank rx.}
#'   \item{Al}{the local low-dimensional structure A, a (jy*ry) by r matrix of rank r.}
#'   \item{Bl}{the local low-dimensional structure B, a (jx*rx) by r matrix of rank r.}
#'   \item{C}{the coefficient matrix C.}
#'   \item{Alstar}{A* in equation (7).}
#'   \item{Blstar}{B* in equation (7).}
#'   \item{Cstar}{\eqn{(U \otimes I_jy)A* B*^T(V \otimes I_jx)^T} in equation (7), the coefficient matrix
#'            used in data generation.}
#'   \item{tseq}{a sequence of the observed time points of y(t) of length nt.}
#'   \item{psi}{the set of basis functions to expand y(t).}
#'   \item{Jpsi}{the correlation matrix of psi, and Jpsihalf is (Jpsi)^(1/2).}
#'   \item{sseq}{a sequence of the observed time points of x(s) of length ns.}
#'   \item{phi}{the set of basis functions to expand x(s).}
#'   \item{Jphi}{the correlation matrix of phi.}
#'   \item{E}{the random error matrix.}
#'   \item{Y}{Y(t), an array of dimension (n, d, nt), with random error.}
#'   \item{X}{X(s), an array of dimension (n, p, ns).}
#'   \item{Ytrue}{Ytrue(t), an array of fimension (n, d, nt), without random error.}
#'   \item{Yest}{response matrix of dimension n by (jy*d), used in estimation.}
#'   \item{Xest}{design matrix of dimension n by (jx*p), used in estimation.}
#' @examples
#' library(NRRR)
#' simDat <- nrrr.sim(n=100,ns=200,nt=200,r=5,rx=3,ry=3,
#'                    jx=15,jy=15,p=10,d=6,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' simDat$Ag
#' @importFrom stats runif rnorm sd
#' @importFrom splines bs
#' @importFrom MASS mvrnorm
#' @export
nrrr.sim <- function(n,ns,nt,r,rx,ry,jx,jy,p,d,s2n,rho_X,rho_E,Sigma=CorrAR){

  #require(MASS)
  #require(splines)

  ### generate X(s) and Y(t)
  # generate uniformly distributed time points for x(s)
  sseq <- sort(stats::runif(ns,0,1))
  # generate basis functions
  phi <- splines::bs(c(0,sseq),df = jx)[-1,]
  # compute Jphi
  Jphi <- matrix(nrow=jx,ncol=jx,0)
  sdiff <- (sseq - c(0,sseq[-ns]))
  for(s in 1:ns) Jphi <- Jphi + phi[s,]%*%t(phi[s,])*sdiff[s]

  # generate X(s)
  Xsigma <- Sigma(p*jx,rho_X)
  Xi <- MASS::mvrnorm(n,rep(0,p*jx),Xsigma)
  X <- array(dim=c(n,p,ns),NA)
  for(s in 1:ns){
    Phis <- kronecker(diag(nrow=p,ncol=p),t(phi[s,]))
    for(i in 1:n){
      X[i,,s] <-  Phis%*%Xi[i,]
    }
  }

  #####Bg (V) matrix, p times rx
  Bg <- qr.Q(qr(matrix(nrow=p,ncol=rx,stats::rnorm(p*rx))))
  #####Ag (U) matrix, d times ry
  Ag <- qr.Q(qr(matrix(nrow=d,ncol=ry,stats::rnorm(d*ry))))
  #####Blstar (B^*) matrix, jx*rx times r
  Blstar <- matrix(nrow=jx*rx, ncol=r, stats::rnorm(jx*rx*r))
  #####Alstar (A^*) matrix, jy*ry times r
  Alstar <- matrix(nrow=jy*ry, ncol=r, stats::rnorm(jy*ry*r))
  #####The coefficient matrix to generate Y(t)
  Cstar <- kronecker(Ag,diag(jy))%*%Alstar%*%t(Blstar)%*%kronecker(t(Bg),diag(jx))


  # Integrated X
  Xint <- matrix(nrow=p*jx,ncol=n,0)
  sdiff <- (sseq - c(0,sseq[-ns]))
  for(s in 1:ns){
    Xint <- Xint + kronecker(diag(nrow=p,ncol=p),phi[s,])%*%t(X[,,s])*sdiff[s]
  }

  # True signal Ytrue(t)
  # generate uniformly distributed time points for y(t)
  tseq <- sort(stats::runif(nt,0,1))
  # generate basis function
  psi <- splines::bs(c(0,tseq),df = jy)[-1,]
  Ytrue <- array(dim=c(n,d,nt),NA)
  Y0 <- Cstar%*%Xint
  for(t in 1:nt){
    Psit <- kronecker(diag(nrow=d,ncol=d),t(psi[t,]))
    Ytrue[,,t] <- t(Psit%*%Y0)
  }


  # The error matrix, n times d*jy
  sigma <- stats::sd(as.vector(Y0))/s2n
  E <- matrix(nrow=n,ncol=d*jy, stats::rnorm(n*d*jy,0,sigma))


  # The response matrix Y(t)
  Y <- array(dim=c(n,d,nt),NA)
  for(t in 1:nt){
    Psit <- kronecker(diag(nrow=d,ncol=d),t(psi[t,]))
    for(i in 1:n){
      Y[i,,t] <- Ytrue[i,,t]+Psit%*%E[i,]
    }
  }


  # compute Jpsi
  Jpsi <- matrix(nrow=jy,ncol=jy,0)
  tdiff <- (tseq - c(0,tseq[-nt]))
  for(t in 1:nt) Jpsi <- Jpsi + psi[t,]%*%t(psi[t,])*tdiff[t]
  eJpsi <- eigen(Jpsi)
  Jpsihalf <- eJpsi$vectors%*%diag(sqrt(eJpsi$values))%*%t(eJpsi$vectors)
  Jpsihalfinv <- eJpsi$vectors%*%diag(1/sqrt(eJpsi$values))%*%t(eJpsi$vectors)


  ### Process the data for estimation
  Xa <- array(dim=c(n,p,jx),NA)
  sdiff <- (sseq - c(0,sseq[-ns]))
  for(i in 1:n){
    for(l in 1:p){
      for(j in 1:jx){
        Xa[i,l,j] <- sum(phi[,j]*X[i,l,]*sdiff)    # integrate over s
      }
    }
  }

  Ya <- array(dim=c(n,d,jy),NA)
  tdiff <- (tseq - c(0,tseq[-nt]))
  psistar <- psi%*%Jpsihalfinv
  for(k in 1:d){
    for(i in 1:n){
      for(j in 1:jy){
        Ya[i,k,j] <- sum(psistar[,j]*Y[i,k,]*tdiff) # integrate over t
      }
    }
  }

  # Y and X matrices for estimation
  Yest <- Ya[,,1]
  for(j in 2:jy) Yest <- cbind(Yest,Ya[,,j])
  Xest <- Xa[,,1]
  for(j in 2:jx) Xest <- cbind(Xest,Xa[,,j])


  # C matrix in estimation
  Alstaradj <-  kronecker(diag(ry),Jpsihalf)%*%Alstar
  alindex <- rep(1:jy,ry)
  Al <- Alstaradj[order(alindex),]
  blindex <- rep(1:jx,rx)
  Bl <- Blstar[order(blindex),]
  C <- kronecker(diag(jx),Bg)%*%Bl%*%t(Al)%*%kronecker(diag(jy),t(Ag))


  list(Ag=Ag,Bg=Bg,Alstar=Alstar,Blstar=Blstar,Cstar=Cstar,
       Al=Al, Bl=Bl, C=C,
       sseq=sseq,psi=psi,Jpsi=Jpsi,Jpsihalf=Jpsihalf,
       tseq=tseq,phi=phi,Jphi=Jphi,
       E=E,
       Y=Y,X=X,Ytrue=Ytrue,
       Yest=Yest,Xest=Xest)

}

