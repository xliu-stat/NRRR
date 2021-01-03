#' @title
#' Select ranks with an information criterion
#'
#' @description
#' This function selects the optimal ranks \code{(r, rx, ry)} using a user-specified
#' information criterion. The degrees of freedom is estimated by the number of
#' free parameters in the model. The blockwise coordinate descent algorithm is used to fit
#' the model with any combinations of \code{(r, rx, ry)}.
#'
#' @usage
#' NestRRR.select(Y, X, Ag0=NULL, Bg0=NULL,
#'                jx, jy, p, d, n, maxiter = 300, conv = 1e-4,
#'                quietly = FALSE, method = c('RRR','RRS')[1],
#'                lambda = 0, ic = c("BIC","BICP","AIC","GCV")[1],
#'                dimred = c(TRUE,TRUE,TRUE), rankfix = NULL)
#'
#' @param Y response matrix of dimension n-by-jy*d.
#' @param X design matrix of dimension n-by-jx*p.
#' @param Ag0 an initial estimator of matrix U. If NULL then generate it
#'            by \code{\link{NestRRRini}}. Default is NULL.
#' @param Bg0 an initial estimator of matrix V, if NULL then generate it
#'            by \code{\link{NestRRRini}}. Default is NULL.
#' @param jx the number of basis functions to expand the functional predictor.
#' @param jy the number of basis functions to expand the functional response.
#' @param p the number of predictors.
#' @param d the number of responses.
#' @param n sample size.
#' @param maxiter the maximum iteration number of the
#'                blockwise coordinate descent algorithm. Default is 300.
#' @param conv the tolerance level used to control the convergence of the
#'             blockwise coordinate descent algorithm. Default is 1e-4.
#' @param quietly a logical value with two options. FALSE (default): show the
#'                rank selection process; TRUE: do not show the process.
#' @param method 'RRR' (default): no additional ridge penalty; 'RRS': add an
#'               additional ridge penalty.
#' @param lambda the tuning parameter to control the amount of ridge
#'               penalization. It is only used when \code{method = 'RRS'}.
#'               Default is 0.
#' @param ic the user-specified information criterion. Four options are available,
#'           including BIC, BICP, AIC, GCV.
#' @param dimred a vector of logical values to decide whether to use the selected information criterion
#'               do rank selection on certain dimensions. TRUE means the rank is selected
#'               by the selected information criterion. If \code{dimred[1]=FALSE}, r is
#'               provided by \code{rankfix} or \eqn{min(jx*p,jy*d, rank(X))};
#'               If \code{dimred[2]=FALSE}, rx equals to p; If \code{dimred[3]=FALSE},
#'               ry equals to d. Default is \code{c(TRUE,TRUE,TRUE)}.
#' @param rankfix a user-provided value of r when \code{dimred[1]=FALSE}. Default is NULL
#'                which leads to \eqn{r=min(jy*d,rank(X))}.
#' @return The function returns a list:
#'   \item{Ag}{the estimated U.}
#'   \item{Bg}{the estimated V.}
#'   \item{Al}{the estimated A.}
#'   \item{Bl}{the estimated B.}
#'   \item{C}{the estimated coefficient matrix C.}
#'   \item{df}{the estimated degrees of freedom of the NRRR model.}
#'   \item{sse}{the sum of squared errors of the selected model.}
#'   \item{ic}{a vector containing values of BIC, BICP, AIC, GCV of the selected model.}
#'   \item{rank}{the estimated r.}
#'   \item{rx}{the estimated rx.}
#'   \item{ry}{the estimated ry.}
#'
#'
#' @details
#' Denote \eqn{\hat C(r, rx, ry)} as the estimator of the coefficient matrix
#' with the rank values fixed at some \eqn{(r, rx, ry)} and write the sum of
#' squared errors as \eqn{SSE(r, rx, ry)=||Y - X\hat C(r, rx, ry)||_F^2}. We
#' define
#' \deqn{BIC(r, rx, ry) = n*d*jy*log(SSE(r, rx, ry)/(n*d*jy)) + log(n*d*jy)*df(r, rx, ry),}
#' where \eqn{df(r, rx, ry)} is the effective degrees of freedom and is estimated
#' by the number of free parameters
#' \deqn{\hat df(r, rx, ry) = rx*(rank(X)/jx - rx) + ry*(d - ry) + (jy*ry + jx*rx - r)*r.}
#' Similarly, we can define the other information criteria. With the defined BIC,
#' a three-dimensional grid search procedure of the rank
#' values is performed, and the best model is chosen as the one with the
#' smallest BIC value. Instead of a nested rank selection method, we apply a
#' one-at-a-time selection approach. We first set \eqn{rx = p, ry = d}, and
#' select the best local rank \eqn{\hat r} among the models with
#' \eqn{1 \le r \le min(rank(X), jy*d)}. We then fix the local rank at
#' \eqn{\hat r} and repeat a similar procedure to determine \eqn{\hat rx}
#' and \eqn{\hat ry}, one at a time. Finally, with fixed \eqn{\hat rx} and \eqn{\hat ry},
#' we refine the estimation of r.
#'
#'
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
# @importFrom rrpack cv.rrr
#' @export
#' @examples
#' library(NRRR)
#' set.seed(3)
#' # Simulation setting 2 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=100,nt=100,r=3,rx=3,ry=3,
#'                    jx=8,jy=8,p=20,d=20,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma="CorrAR")
#' fit_init <- with(simDat, NestRRR.select(Yest,Xest,Ag0=NULL,Bg0=NULL,
#'                               jx=8,jy=8,p=20,d=20,n=100,
#'                               maxiter=300,conv=1e-4,quietly=FALSE,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),rankfix=NULL))
NestRRR.select <- function(Y,X,Ag0=NULL,Bg0=NULL,
                           jx,jy,p,d,n,maxiter=300,conv=1e-4,
                           quietly=FALSE,method=c('RRR','RRS')[1],
                           lambda=0,ic=c("BIC","BICP","AIC","GCV")[1],
                           dimred = c(TRUE,TRUE,TRUE),rankfix=NULL
){
  #require(rrpack)
  if (method == "RRS" & lambda == 0) stop("A positive tuning parameter should be provided when 'RRS' is used.")
  if (!(ic %in% c("BIC","BICP","AIC","GCV"))) stop("A valid information criterion name should be provided.")

  # compute rank(X)
  xr <- sum(svd(X)$d>1e-2)

  # initialize r
  if(dimred[1]){
    fitRRR <- cv.rrr(Y,X,nfold=10)
    rest <- fitRRR$rank
    # If zero fit
    if(rest==0){
      fitRRR <- RRR(Y,X,nrank=1)
      rest <- fitRRR$rank
    }
    if(!quietly) {
      cat("Initial r   = ",rest, "\n",sep="")
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
