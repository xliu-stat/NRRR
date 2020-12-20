#' @title
#' Generate initial estimators
#'
#' @description
#' This function provides the initial estimators of U and V to initialize the
#' blockwise coordinate descent algorithm.
#'
#' @usage
#' NestRRRini(Y, X, r, rx, ry, jx, jy, p, d, n)
#'
#'
#' @param Y the response matrix of dimension n-by-jy*d.
#' @param X the design matrix of dimension n-by-jx*p.
#' @param r rank of the local reduced-rank structure.
#' @param rx the number of latent predictors.
#' @param ry the number of latent responses.
#' @param jx the number of basis functions to expand functional predictor.
#' @param jy the number of basis functions to expand functional response.
#' @param p the number of predictors.
#' @param d the number of responses.
#' @param n the sample size.
#'
#'
#' @return The function returns a list:
#'   \item{Ag}{the estimated U.}
#'   \item{Bg}{the estimated V.}
#'
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#'
#' @examples
#' library(NRRR)
#' simDat <- nrrr.sim(n=100,ns=200,nt=200,r=5,rx=3,ry=3,
#'                    jx=15,jy=15,p=10,d=6,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_init <- with(simDat, NestRRRini(Y=Yest,X=Xest,r=5,
#'                  rx=3,ry=3,jx=15,jy=15,p=10,d=6,n=100))
#' fit_init$Ag
#' @export
NestRRRini <- function(Y,X,r,rx,ry,jx,jy,p,d,n){
  #require(rrpack)
  # ignore the global structure to estimate C from RRR
  fit_RRR <- RRR(Y,X,nrank=r)
  Crr <- fit_RRR$C
  # Compute Bghat, i.e, V
  if(p==rx){
    Bghat <- diag(nrow=rx,ncol=rx)
  }else{
    Cbg <- matrix(nrow=p,ncol=d*jy*jx,NA)
    for(j in 1:jx){
      a1 <- d*jy*(j-1) + 1
      b1 <- d*jy*j

      a2 <- p*(j-1) + 1
      b2 <- p*j
      Cbg[,a1:b1] <- Crr[a2:b2,]
    }
    Bghat <- as.matrix(svd(Cbg, nu=rx,nv=rx)$u)
  }


  #Compute Aghat, i.e., U
  if(d==ry){
    Aghat <- diag(nrow=ry,ncol=ry)
  }else{
    Cag <- matrix(nrow=d,ncol=p*jy*jx,NA)
    for(j in 1:jy){
      a1 <- p*jx*(j-1) + 1
      b1 <- p*jx*j

      a2 <- d*(j-1) + 1
      b2 <- d*j

      Cag[,a1:b1] <- t(Crr[,a2:b2])
    }
    Aghat <- as.matrix(svd(Cag, nu=ry,nv=ry)$u)
  }


  return(list(Ag=Aghat,Bg=Bghat,C=Crr))
}



