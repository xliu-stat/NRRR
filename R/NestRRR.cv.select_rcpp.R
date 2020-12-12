#' NRRR with cross validation to select ranks (Rcpp)
#'
#' This function applies Cross-Validation
#'  to select the optimal (r, rx, ry) with
#'  the blockwise coordinate descent algorithm.
#'
#' @param Y Response matrix.
#' @param X Design matrix.
#' @param nfold Number of folds, default is 10.
#' @param norder Assign samples to multiple folds for cross validation.
#' @param jx Number of basis functions to expand x(s).
#' @param jy Number of basis functions to expand y(t).
#' @param p Number of predictors.
#' @param d Number of responses.
#' @param n Sample size.
#' @param maxiter Maximum iteration number, Default 300.
#' @param conv Tolerance level to control convergence, default 1e-4.
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
#'   \item{rank}{a scalar, the estimated r.}
#'   \item{rx}{a scalar, the estimated rx.}
#'   \item{ry}{a scalar, the estimated ry.}
#'   \item{rxErrmat}{a matrix displays errors in the path of selecting rx.
#'                   For each element, 0 indicates no error.}
#'   \item{ryErrmat}{a matrix displays errors in the path of selecting ry.
#'                   For each element, 0 indicates no error.}
#'   \item{rErrmat}{a matrix displays errors in the path of selecting r.
#'                   For each element, 0 indicates no error.}
#' @references Liu, X., Ma, S., & Chen, K. (2020).
#' Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#' @export
#' @examples
#' library(NRRR)
#' set.seed(1)
#' # Simulation setting 1 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=60,nt=60,r=5,rx=3,ry=3,
#'                    jx=8,jy=8,p=10,d=10,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma=CorrAR)
#' fit_init <- with(simDat, NestRRR.cv.select_rcpp(Yest,Xest,nfold=10,norder=NULL,
#'                               jx=8,jy=8,p=10,d=10,n=100,
#'                               maxiter=300,conv=1e-4,
#'                               method=c('RRR','RRS')[1],lambda=0,
#'                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),
#'                               rankfix=NULL,xrankfix=NULL,yrankfix=NULL))
NestRRR.cv.select_rcpp <- function(Y,X,nfold=10,norder=NULL,jx,jy,p,d,n,
                              maxiter=300,conv=1e-4,method=c('RRR','RRS')[1],lambda=0,
                              ic=c("BIC","BICP","AIC","GCV")[1],
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


  if (ic == 'BIC'){
    ic_num <- 0
  } else if (ic == 'BICP'){
    ic_num <- 1
  } else if (ic == 'AIC'){
    ic_num <- 2
  } else {
    ic_num <- 3
  }


  method <- ifelse(method == 'RRR', 1, 2)

  dimred1 <- ifelse(dimred[1],1,0)
  dimred2 <- ifelse(dimred[2],1,0)
  dimred3 <- ifelse(dimred[3],1,0)

  xrankfix <- ifelse(is.null(xrankfix),0,xrankfix)
  yrankfix <- ifelse(is.null(yrankfix),0,yrankfix)

  norder <- norder - 1

  fit <- nrrr_cv_my(Y, X, norder, nfold, xr, rfit, xrankfix, yrankfix,
                    jx, jy, p, d, n, ic_num, maxiter, method,
                    dimred1, dimred2, dimred3, conv, lambda)

  return(list(Ag=fit$Ag,Bg=fit$Bg,Al=fit$Al,Bl=fit$Bl,C=fit$C,df=fit$df,
              sse=fit$sse,ic=fit$ic,obj=fit$obj,rx_path=fit$rx_path,ry_path=fit$ry_path,
              r_path=fit$r_path,
              rank=fit$rank,rx=fit$rx,ry=fit$ry,
              rxErrmat=fit$rxErrmat,ryErrmat=fit$ryErrmat,rErrmat=fit$rErrmat))
}
