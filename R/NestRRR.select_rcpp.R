#' Select ranks with information criterion (Rcpp)
#'
#' This function select the optimal (r, rx, ry) with a specified information criterion
#' based on the proposed blockwise coordinate descent algorithm.
#'
#' @param Y Response matrix with n rows and jy*d columns.
#' @param X Design matrix with n rows and jx*p columns.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param n Sample size.
#' @param maxiter Maximum iteration number of the blockwise coordinate descent algorithm, default 300.
#' @param conv Tolerance level to control convergence, default 1e-4.
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
#'           which leads to min(jx*p,jy*d,rank(X)).
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
#'   \item{rxErrseq}{a vector displays errors in the path of selecting rx.
#'                   For each element, 0 indicates no error.}
#'   \item{ryErrseq}{a vector displays errors in the path of selecting ry.
#'                   For each element, 0 indicates no error.}
#'   \item{rErrseq}{a vector displays errors in the path of selecting r.
#'                   For each element, 0 indicates no error.}

#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#' @export
#' @examples
#' library(NRRR)
#' set.seed(3)
#' # Simulation setting 2 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=100,nt=100,r=3,rx=3,ry=3,
#'                    jx=8,jy=8,p=20,d=20,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_init <- with(simDat, NestRRR.select_rcpp(Yest,Xest,
#'                               jx=8,jy=8,p=20,d=20,n=100,
#'                               maxiter=300,conv=1e-4,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),rankfix=NULL))
NestRRR.select_rcpp <- function(Y,X,jx,jy,p,d,n,maxiter=300,
                                conv=1e-4,method=c('RRR','RRS')[1],
                                lambda=0,ic=c("BIC","BICP","AIC","GCV")[1],
                                ##if do not do rank reduction, then use rankfix
                                dimred = c(TRUE,TRUE,TRUE),rankfix=NULL
){
  #require(rrpack)
  # compute rank(X)
  xr <- sum(svd(X)$d>1e-2)

  method <- ifelse(method == 'RRR', 1, 2)

  dimred1 <- ifelse(dimred[1],1,0)
  dimred2 <- ifelse(dimred[2],1,0)
  dimred3 <- ifelse(dimred[3],1,0)


  if (ic == 'BIC'){
    ic_num <- 0
  } else if (ic == 'BICP'){
    ic_num <- 1
  } else if (ic == 'AIC'){
    ic_num <- 2
  } else {
    ic_num <- 3
  }

  # initialize r
  if(dimred[1]){
    fitRRR <- rrpack::cv.rrr(Y,X,nfold=10)
    rest <- fitRRR$rank
    # if(!quietly) {
    #   cat("Initial r   = ",rest, "\n",sep="")
    # }
    # If zero fit
    if(rest==0){
      fitRRR <- RRR(Y,X,nrank=1)
      rest <- fitRRR$rank
    }
  }else{
    rest <- ifelse(is.null(rankfix),min(ncol(Y),ncol(X),xr),rankfix)
  }
  rfit <- rest

  fit <- nrrr_select_my(Y, X, xr, rfit, jx, jy, p, d, n, ic_num, maxiter, conv,
                        method, lambda, dimred1, dimred2, dimred3)

  return(list(Ag=fit$Ag,Bg=fit$Bg,Al=fit$Al,Bl=fit$Bl,C=fit$C,df=fit$df,
              sse=fit$sse,ic=fit$ic,obj=fit$obj,
              rank=fit$rank,rx=fit$rx,ry=fit$ry,
              rxErrseq=fit$rxErrseq,ryErrseq=fit$ryErrseq,rErrseq=fit$rErrseq))
}
