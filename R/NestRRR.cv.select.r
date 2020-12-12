#' NRRR with cross validation to select ranks
#'
#' This function applies Cross-Validation
#'  to select the optimal (r, rx, ry) with
#'  the blockwise coordinate descent algorithm.
#'
#' @param Y Response matrix.
#' @param X Design matrix.
#' @param nfold Number of folds, default is 10.
#' @param norder Assign samples to multiple folds for cross validation.
#' @param Ag0 Initial estimator of U, if NULL then generate it by
#'           function \code{NestRRRini()}.
#' @param Bg0 Initial estimator of V, if NULL then generate it by
#'           function \code{NestRRRini()}.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param n Sample size.
#' @param maxiter Maximum iteration number, Default 300.
#' @param conv Tolerance level to control convergence, default 1e-4.
#' @param quietly FALSE (default) show the rank selection process; TRUE not show.
#' @param method RRR (default) no ridge penalty; RRS use ridge penalty.
#' @param lambda Tuning parameter for the ridge penalty, only used when
#'           method='RRS', default=0.
#' @param ic Information criterion, selected from BIC, BICP, AIC, GCV.
#' @param dimred A vector to decide whether do dimension reduction on certain
#'           dimensions, with c(TRUE,TRUE,TRUE) corresponding to c(r,rx,ry),
#'           TRUE (default) do dimension reduction; FALSE do not.
#'           If dimred[1]=FALSE, r is provided by rankfix or min(jx*p,jy*d,rank(X));
#'           If dimred[2]=FALSE, rx is provided by xrankfix or p;
#'           If dimred[3]=FALSE, ry is provided by yrankfix or d.
#' @param rankfix Provide a value for r when dimred[1]=FALSE, default is NULL
#'           which leads to min(jx*p,jy*d,rank(X)).
#' @param xrankfix Provide a value for rx when dimred[2]=FALSE, default is NULL
#'           which leads to p.
#' @param yrankfix Provide a value for ry when dimred[3]=FALSE, default is NULL
#'           which leads to d.
#' @return The returned items are
#'   \item{Ag}{estimated matrix U.}
#'   \item{Bg}{estimated matrix V.}
#'   \item{Al}{estimated matrix A.}
#'   \item{Bl}{estimated matrix B.}
#'   \item{C}{estimated matrix C.}
#'   \item{df}{a scalar, estimated degrees of freedom of the selected model.}
#'   \item{sse}{a scalar, the sum of squared errors of the selected model.}
#'   \item{ic}{a vector contains values of BIC,BICP,AIC,GCV of the selected model.}
#'   \item{obj}{a vector contains all objective function (sse) values along
#'           iterations of the selected model.}
#'   \item{rx_path}{a matrix displays the path of selecting rx with cross validation.}
#'   \item{ry_path}{a matrix displays the path of selecting ry with cross validation.}
#'   \item{r_path}{a matrix displays the path of selecting r with cross validation.}
#'   \item{r_fitseq}{a sequence of possible r values.}
#'   \item{iter}{the number of iterations needed to converge in the selected model.}
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
#' set.seed(1)
#' # Simulation setting 1 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=60,nt=60,r=5,rx=3,ry=3,
#'                    jx=8,jy=8,p=10,d=10,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_init <- with(simDat, NestRRR.cv.select(Yest,Xest,nfold=10,norder=NULL,
#'                               Ag0=NULL,Bg0=NULL,
#'                               jx=8,jy=8,p=10,d=10,n=100,
#'                               maxiter=300,conv=1e-4,quietly=FALSE,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),
#'                               rankfix=NULL,xrankfix=NULL,yrankfix=NULL))
NestRRR.cv.select <- function(Y,X,nfold=10,norder=NULL,Ag0=NULL,Bg0=NULL,jx,jy,p,d,n,
                              maxiter=300,conv=1e-4,quietly=FALSE,
                              method=c('RRR','RRS')[1],lambda=0,
                              ic=c("BIC","BICP","AIC","GCV")[1],
                              ##if do not do rank reduction, then use rankfix
                              dimred = c(TRUE,TRUE,TRUE),
                              rankfix=NULL,xrankfix=NULL,yrankfix=NULL
){
  #require(rrpack)
  xr <- sum(svd(X)$d>1e-2)
  if (is.null(norder))
    norder <- sample(seq_len(n),n)

  # initialize r
  if(dimred[1]){
    fitRRR <- rrpack::cv.rrr(Y,X,nfold=10,norder = norder)
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

  # Select rx by cross validation
  rx_path <- matrix(ncol = nfold, nrow = p, NA)
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

  if(!quietly) {
    cat("Selected rx = ",rxest, "\n",sep="")
  }

  # Select ry by cross validation
  ry_path <- matrix(ncol = nfold, nrow = d, NA)
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
