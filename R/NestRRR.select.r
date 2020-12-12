#' Select ranks with information criterion
#'
#' This function select the optimal (r, rx, ry) with a specified information criterion
#' based on the proposed blockwise coordinate descent algorithm.
#'
#' @param Y Response matrix with n rows and jy*d columns.
#' @param X Design matrix with n rows and jx*p columns.
#' @param Ag0 Initial estimator of U, if NULL then generate it by
#'           function \code{NestRRRini()}.
#' @param Bg0 Initial estimator of V, if NULL then generate it by
#'           function \code{NestRRRini()}.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param n Sample size.
#' @param maxiter Maximum iteration number of the blockwise coordinate descent algorithm, default 300.
#' @param conv Tolerance level to control convergence, default 1e-4.
#' @param quietly FALSE (default) show the rank selection process, if TRUE then not.
#' @param method RRR (default) no ridge penalty; RRS use ridge penalty.
#' @param lambda Tuning parameter for the ridge penalty, only used when method='RRS', default=0.
#' @param ic Information criterion, selected from BIC, BICP, AIC, GCV.
#' @param dimred Vector to decide whether do dimension reduction on certain
#'           dimensions, with c(TRUE,TRUE,TRUE) corresponding to c(r,rx,ry),
#'           TRUE (default) do dimension reduction; FALSE do not.
#'           If dimred[1]=FALSE, r is provided by rankfix or min(jx*p,jy*d, rank(X));
#'           If dimred[2]=FALSE, rx equals to p;
#'           If dimred[3]=FALSE, ry equals to d.
#' @param rankfix Provide a value for r when dimred[1]=FALSE, default is NULL
#'           which leads to min(jx*p,jy*d,rank(X))
#' @return The returned results containing
#'   \item{Ag}{the global low-dimensional structure U, a d by ry matrix of rank ry.}
#'   \item{Bg}{the global low-dimensional structure V, a p by rx matrix of rank rx.}
#'   \item{Al}{the local low-dimensional structure A.}
#'   \item{Bl}{the local low-dimensional structure B.}
#'   \item{C}{the NRRR estimator of the coefficient matrix C.}
#'   \item{df}{a scalar, the estimated degrees of freedom of the NRRR model.}
#'   \item{sse}{a scalar, the sum of squared errors of the selected model.}
#'   \item{ic}{a vector contains values of BIC,BICP,AIC,GCV of the selected model.}
#'   \item{obj}{ a vector contains all objective function (sse) values along iterations of the selected model.}
#'   \item{rank}{a scalar, the estimated r.}
#'   \item{rx}{a scalar, the estimated rx.}
#'   \item{ry}{a scalar, the estimated ry.}
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#' @importFrom rrpack cv.rrr
#' @export
#' @examples
#' library(NRRR)
#' set.seed(3)
#' # Simulation setting 2 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=100,nt=100,r=3,rx=3,ry=3,
#'                    jx=8,jy=8,p=20,d=20,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_init <- with(simDat, NestRRR.select(Yest,Xest,Ag0=NULL,Bg0=NULL,
#'                               jx=8,jy=8,p=20,d=20,n=100,
#'                               maxiter=300,conv=1e-4,quietly=FALSE,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),rankfix=NULL))
NestRRR.select <- function(Y,X,Ag0=NULL,Bg0=NULL,jx,jy,p,d,n,maxiter=300,
                           conv=1e-4,quietly=FALSE,method=c('RRR','RRS')[1],
                           lambda=0,ic=c("BIC","BICP","AIC","GCV")[1],
                           ##if do not do rank reduction, then use rankfix
                           dimred = c(TRUE,TRUE,TRUE),rankfix=NULL
){
  #require(rrpack)
  # compute rank(X)
  xr <- sum(svd(X)$d>1e-2)

  # initialize r
  if(dimred[1]){
    fitRRR <- rrpack::cv.rrr(Y,X,nfold=10)
    rest <- fitRRR$rank
    if(!quietly) {
      cat("Initial r   = ",rest, "\n",sep="")
    }
    # If zero fit
    if(rest==0){
      fitRRR <- RRR(Y,X,nrank=1)
      rest <- fitRRR$rank
    }
  }else{
    rest <- ifelse(is.null(rankfix),min(ncol(Y),ncol(X),xr),rankfix)
  }
  rfit <- rest


  # Select rx
  if(dimred[2]){
    rxfitseq <- 1:p
    ryfitseq <- rep(d,length(rxfitseq))

    icseq <- vector()
    for(i in 1:length(rxfitseq)){
      rxfit <- rxfitseq[i]
      ryfit <- ryfitseq[i]
      fit1 <- NestRRR(Y,X,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,p,
                      d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                      method=method,lambda=lambda)
      icseq[i] <- switch(ic,'BIC'=fit1$ic[1],'BICP'=fit1$ic[2],
                         'AIC'=fit1$ic[3],'GCV'=fit1$ic[4])
    }
    rxest <- rxfitseq[which.min(icseq)]
  }else{
    rxest <- p
  }
  rxfit <- rxest
  if(!quietly) {
    cat("Selected rx = ",rxest, "\n",sep="")
  }


  # Select ry
  if(dimred[3]){
    ryfitseq <- 1:d
    rxfitseq <- rep(rxest,length(ryfitseq))

    icseq <- vector()
    for(i in 1:length(rxfitseq)){
      rxfit <- rxfitseq[i]
      ryfit <- ryfitseq[i]
      fit1 <- NestRRR(Y,X,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,
                      p,d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                      method=method,lambda=lambda)
      icseq[i] <- switch(ic,'BIC'=fit1$ic[1],'BICP'=fit1$ic[2],
                         'AIC'=fit1$ic[3],'GCV'=fit1$ic[4])
    }
    ryest <- ryfitseq[which.min(icseq)]
  }else{
    ryest <- d
  }
  ryfit <- ryest
  if(!quietly) {
    cat("Selected ry = ",ryest, "\n",sep="")
  }

  # refine rank selection
  if(dimred[1]){
    rxfit <- rxest
    ryfit <- ryest
    icseq <- vector()
    # range of the possible rank
    rfitseq <- max(1,rest-5):min(rest+5,min(n,p*jx,d*jy))

    rfit <- rfitseq[1]
    fit <- NestRRR(Y,X,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,
                   p,d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                   method=method,lambda=lambda)
    icseq[1] <- switch(ic,'BIC'=fit$ic[1],'BICP'=fit$ic[2],
                       'AIC'=fit$ic[3],'GCV'=fit$ic[4])
    icmin <- icseq[1]
    for(i in 2:length(rfitseq)){
      rfit <- rfitseq[i]
      fit1 <- NestRRR(Y,X,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,
                      p,d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                      method=method,lambda=lambda)
      icseq[i] <- switch(ic,'BIC'=fit1$ic[1],'BICP'=fit1$ic[2],
                         'AIC'=fit1$ic[3],'GCV'=fit1$ic[4])
      if(icseq[i]<icmin){
        fit <- fit1
        icmin <- icseq[i]
      }
    }
    rest <- rfitseq[which.min(icseq)]
  }else{
    fit <- NestRRR(Y,X,NULL,NULL,rini=rest,rest,rxfit,ryfit,jx,jy,
                   p,d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                   method=method,lambda=lambda)
  }

  if(!quietly) {
    cat("Selected r  = ",rest,"\n", sep="")
  }

  return(list(Ag=fit$Ag,Bg=fit$Bg,Al=fit$Al,Bl=fit$Bl,C=fit$C,df=fit$df,
              sse=fit$sse,ic=fit$ic,obj=fit$obj,
              rank=rest,rx=rxest,ry=ryest))

}
