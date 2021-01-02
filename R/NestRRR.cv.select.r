#' @title
#' Select ranks with cross validation
#'
#' @description
#' This function selects the optimal ranks \code{(r, rx, ry)} using a cross
#' validation procedure. The blockwise coordinate descent algorithm is used to fit
#' the model with any combinations of \code{(r, rx, ry)}.
#'
#' @usage
#' NestRRR.cv.select(Y, X, nfold = 10, norder = NULL, Ag0 = NULL, Bg0 = NULL,
#'                   jx, jy, p, d, n,
#'                   maxiter = 300, conv = 1e-4, quietly = FALSE,
#'                   method = c('RRR','RRS')[1], lambda=0,
#                   ic=c("BIC","BICP","AIC","GCV")[1],
#'                   dimred = c(TRUE,TRUE,TRUE),
#'                   rankfix = NULL, xrankfix = NULL, yrankfix = NULL)
#'
#' @param Y response matrix of dimension n-by-jy*d.
#' @param X design matrix of dimension n-by-jx*p.
#' @param nfold the number of folds used in cross validation. Default is 10.
#' @param norder a vector of length n that assigns samples to multiple folds for cross validation.
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
# @param ic the user-specified information criterion. Four options are available,
#'           including BIC, BICP, AIC, GCV.
#' @param dimred a vector of logical values to decide whether to use cross validation
#'               do rank selection on certain dimensions. TRUE (default): yes; FALSE: no.
#'               If \code{dimred[1]=FALSE}, r is provided by \code{rankfix}
#'               or \eqn{min(jy*d, rank(X))};
#'               If \code{dimred[2]=FALSE}, rx equals to \code{xrankfix} or p; If \code{dimred[3]=FALSE},
#'               ry equals to \code{yrankfix} or d. Default is \code{c(TRUE,TRUE,TRUE)}.
#' @param rankfix a user-provided value of r when \code{dimred[1]=FALSE}. Default is NULL
#'                which leads to \eqn{r=min(jy*d,rank(X))}.
#' @param xrankfix a user-provided value of rx when \code{dimred[2]=FALSE}. Default is NULL
#'                 which leads to \code{rx=p}.
#' @param yrankfix a user-provided value of ry when \code{dimred[3]=FALSE}. Default is NULL
#'                 which leads to \code{ry=d}.
#'
#'
#' @return The function returns a list:
#'   \item{Ag}{the estimated U.}
#'   \item{Bg}{the estimated V.}
#'   \item{Al}{the estimated A.}
#'   \item{Bl}{the estimated B.}
#'   \item{C}{the estimated coefficient matrix C.}
#'   \item{df}{the estimated degrees of freedom of the selected model.}
#'   \item{sse}{the sum of squared errors of the selected model.}
#'   \item{ic}{a vector containing values of BIC, BICP, AIC, GCV of the selected model.}
#'   \item{rx_path}{a matrix displays the path of selecting rx with cross validation.}
#'   \item{ry_path}{a matrix displays the path of selecting ry with cross validation.}
#'   \item{r_path}{a matrix displays the path of selecting r with cross validation.}
#   \item{r_fitseq}{a sequence of possible r values.}
#'   \item{iter}{the number of iterations needed to converge in the selected model.}
#'   \item{rank}{the estimated r.}
#'   \item{rx}{the estimated rx.}
#'   \item{ry}{the estimated ry.}
#'
#' @details
#' A three-dimensional grid search procedure of the rank
#' values is performed, and the best model is chosen as the one with the
#' smallest prediction error. Instead of a nested rank selection method, we apply a
#' one-at-a-time selection approach. We first set \eqn{rx = p, ry = d}, and
#' select the best local rank \eqn{\hat r} among the models with
#' \eqn{1 \le r \le min(rank(X), jy*d)}. We then fix the local rank at
#' \eqn{\hat r} and repeat a similar procedure to determine \eqn{\hat rx}
#' and \eqn{\hat ry}, one at a time. Finally, with fixed \eqn{\hat rx} and \eqn{\hat ry},
#' we refine the estimation of r.
#'
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#'
#' @importFrom rrpack cv.rrr
#' @export
#' @examples
#' library(NRRR)
#' set.seed(1)
#' # Simulation setting 1 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=60,nt=60,r=5,rx=3,ry=3,
#'                    jx=8,jy=8,p=10,d=10,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma="CorrAR")
#' fit_init <- with(simDat, NestRRR.cv.select(Yest,Xest,nfold=10,norder=NULL,
#'                               Ag0=NULL,Bg0=NULL,
#'                               jx=8,jy=8,p=10,d=10,n=100,
#'                               maxiter=300,conv=1e-4,quietly=FALSE,
#'                               method=c('RRR','RRS')[1],lambda=0,
#                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),
#'                               rankfix=NULL,xrankfix=NULL,yrankfix=NULL))
NestRRR.cv.select <- function(Y,X,nfold=10,norder=NULL,Ag0=NULL,Bg0=NULL,jx,jy,p,d,n,
                              maxiter=300,conv=1e-4,quietly=FALSE,
                              method=c('RRR','RRS')[1],lambda=0,
#                              ic=c("BIC","BICP","AIC","GCV")[1],
                              ##if do not do rank reduction, then use rankfix
                              dimred = c(TRUE,TRUE,TRUE),
                              rankfix=NULL,xrankfix=NULL,yrankfix=NULL
){
  #require(rrpack)
  if (method == "RRS" & lambda == 0) stop("A positive tuning parameter should be provided when 'RRS' is used.")

  xr <- sum(svd(X)$d>1e-2)
  if (is.null(norder))
    norder <- sample(seq_len(n),n)

  # initialize r
  if(dimred[1]){
    fitRRR <- rrpack::cv.rrr(Y,X,nfold=10,norder = norder)
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

  # Select rx by cross validation
  rx_path <- matrix(ncol = nfold, nrow = p, NA)
  if (p == 1) {
    rxest <- p
  } else {
    if(dimred[2]){
      ndel <- round(n/nfold)
      for (f in seq_len(nfold)){
        if (f != nfold) {
          iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
        }
        else {
          iddel <- norder[(1 + ndel * (f - 1)):n]
        }
        ndel <- length(iddel)
        nf <- n - ndel
        idkeep <- (seq_len(n))[-iddel]
        Xf <- X[-iddel, ]
        Xfdel <- X[iddel, ]
        Yf <- Y[-iddel, ]
        Yfdel <- Y[iddel, ]
        for (i in seq_len(p)){
          rxfit <- i
          ryfit <- d
          fit1 <- NestRRR(Yf,Xf,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,p,
                          d,n=nf,maxiter=maxiter,conv=conv,quietly=TRUE,
                          method=method,lambda=lambda)
          rx_path[i,f] <- sum((Yfdel-Xfdel%*%fit1$C)^2)
        }
      }
      index <- order(colSums(rx_path))
      crerr <- rowSums(rx_path[, index])/length(index) * nfold
      rxest <- which.min(crerr)
    } else {
      rxest <- ifelse(is.null(xrankfix),p,xrankfix)
    }
  }

  if(!quietly) {
    cat("Selected rx = ",rxest, "\n",sep="")
  }

  # Select ry by cross validation
  ry_path <- matrix(ncol = nfold, nrow = d, NA)
  if (d == 1) {
    ryest <- d
  } else {
    if(dimred[3]){
      ndel <- round(n/nfold)
      for (f in seq_len(nfold)){
        if (f != nfold) {
          iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
        }
        else {
          iddel <- norder[(1 + ndel * (f - 1)):n]
        }
        ndel <- length(iddel)
        nf <- n - ndel
        idkeep <- (seq_len(n))[-iddel]
        Xf <- X[-iddel, ]
        Xfdel <- X[iddel, ]
        Yf <- Y[-iddel, ]
        Yfdel <- Y[iddel, ]
        for (i in seq_len(d)){
          rxfit <- rxest
          ryfit <- i
          fit1 <- NestRRR(Yf,Xf,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,p,
                          d,n=nf,maxiter=maxiter,conv=conv,quietly=TRUE,
                          method=method,lambda=lambda)
          ry_path[i,f] <- sum((Yfdel-Xfdel%*%fit1$C)^2)
        }
      }
      index <- order(colSums(ry_path))
      crerr <- rowSums(ry_path[, index])/length(index) * nfold
      ryest <- which.min(crerr)
    } else {
      ryest <- ifelse(is.null(yrankfix),d,yrankfix)
    }
  }

  if(!quietly) {
    cat("Selected ry = ",ryest, "\n",sep="")
  }


  # refine rank selection by cross validation
  rfitseq <- max(1,rest-5):min(rest+5,min(n,p*jx,d*jy))
  r_path <- matrix(ncol = nfold, nrow = length(rfitseq), NA)
  if(dimred[1]){
    rxfit <- rxest
    ryfit <- ryest

    ndel <- round(n/nfold)
    for (f in seq_len(nfold)){
      if (f != nfold) {
        iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
      }
      else {
        iddel <- norder[(1 + ndel * (f - 1)):n]
      }
      ndel <- length(iddel)
      nf <- n - ndel
      idkeep <- (seq_len(n))[-iddel]
      Xf <- X[-iddel, ]
      Xfdel <- X[iddel, ]
      Yf <- Y[-iddel, ]
      Yfdel <- Y[iddel, ]
      for (i in 1:length(rfitseq)){
        rfit <- rfitseq[i]
        fit1 <- NestRRR(Yf,Xf,NULL,NULL,rini=rfit,rfit,rxfit,ryfit,jx,jy,p,
                        d,n=nf,maxiter=maxiter,conv=conv,quietly=TRUE,
                        method=method,lambda=lambda)
        r_path[i,f] <- sum((Yfdel-Xfdel%*%fit1$C)^2)
      }
    }
    index <- order(colSums(r_path))
    crerr <- rowSums(r_path[, index])/length(index) * nfold
    rest <- rfitseq[which.min(crerr)]
    fit <- NestRRR(Y,X,NULL,NULL,rini=rest,rest,rxfit,ryfit,jx,jy,p,
                   d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                   method=method,lambda=lambda)

  } else {
    rxfit <- rxest
    ryfit <- ryest
    fit <- NestRRR(Y,X,NULL,NULL,rini=rest,rest,rxfit,ryfit,jx,jy,
                   p,d,n,maxiter=maxiter,conv=conv,quietly=TRUE,
                   method=method,lambda=lambda)

  }

  if(!quietly) {
    cat("Selected r  = ",rest,"\n", sep="")
  }
  return(list(Ag=fit$Ag,Bg=fit$Bg,Al=fit$Al,Bl=fit$Bl,C=fit$C,df=fit$df,
              sse=fit$sse,ic=fit$ic,obj=fit$obj,rx_path=rx_path,ry_path=ry_path,
              r_path=r_path,rfitseq=rfitseq,
              iter=fit$iter,
              rank=rest,rx=rxest,ry=ryest))
}
