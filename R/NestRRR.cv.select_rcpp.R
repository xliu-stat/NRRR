#' @title
#' Select ranks with cross validation (Rcpp)
#'
#' @description
#' This function selects the optimal ranks \code{(r, rx, ry)} using a cross
#' validation procedure. The blockwise coordinate descent algorithm is used to fit
#' the model with any combinations of \code{(r, rx, ry)}.
#'
#' @usage
#' NestRRR.cv.select_rcpp(Y, X, nfold = 10, norder = NULL,
#'                        jx, jy, p, d, n,
#'                        maxiter = 300, conv = 1e-4,
#'                        method = c('RRR','RRS')[1], lambda=0,
#                         ic=c("BIC","BICP","AIC","GCV")[1],
#'                        dimred = c(TRUE,TRUE,TRUE),
#'                        rankfix = NULL, xrankfix = NULL, yrankfix = NULL)
#'
#' @param Y response matrix of dimension n-by-jy*d.
#' @param X design matrix of dimension n-by-jx*p.
#' @param nfold the number of folds used in cross validation. Default is 10.
#' @param norder a vector of length n that assigns samples to multiple folds for cross validation.
#' @param jx the number of basis functions to expand the functional predictor.
#' @param jy the number of basis functions to expand the functional response.
#' @param p the number of predictors.
#' @param d the number of responses.
#' @param n sample size.
#' @param maxiter the maximum iteration number of the
#'                blockwise coordinate descent algorithm. Default is 300.
#' @param conv the tolerance level used to control the convergence of the
#'             blockwise coordinate descent algorithm. Default is 1e-4.
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
#'   \item{rank}{the estimated r.}
#'   \item{rx}{the estimated rx.}
#'   \item{ry}{the estimated ry.}
#   \item{obj}{a vector contains all objective function (sse) values along
#           iterations of the selected model.}
#'   \item{rxErrmat}{a matrix of error flags for the path of selecting rx.
#'                   For each element, 0 indicates no error.}
#'   \item{ryErrmat}{a matrix of error flags for the path of selecting ry.
#'                   For each element, 0 indicates no error.}
#'   \item{rErrmat}{a matrix of error flags for the path of selecting r.
#'                   For each element, 0 indicates no error.}
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
#' @export
#' @examples
#' library(NRRR)
#' set.seed(1)
#' # Simulation setting 1 in NRRR paper
#' simDat <- nrrr.sim(n=100,ns=60,nt=60,r=5,rx=3,ry=3,
#'                    jx=8,jy=8,p=10,d=10,s2n=1,rho_X=0.5,
#'                    rho_E=0,Sigma="CorrAR")
#' fit_init <- with(simDat, NestRRR.cv.select_rcpp(Yest,Xest,nfold=10,norder=NULL,
#'                               jx=8,jy=8,p=10,d=10,n=100,
#'                               maxiter=300,conv=1e-4,
#'                               method=c('RRR','RRS')[1],lambda=0,
#                               ic=c("BIC","BICP","AIC","GCV")[1],
#'                               dimred = c(TRUE,TRUE,TRUE),
#'                               rankfix=NULL,xrankfix=NULL,yrankfix=NULL))
NestRRR.cv.select_rcpp <- function(Y,X,nfold=10,norder=NULL,jx,jy,p,d,n,
                              maxiter=300,conv=1e-4,method=c('RRR','RRS')[1],lambda=0,
                              #ic=c("BIC","BICP","AIC","GCV")[1],
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


  method <- ifelse(method == 'RRR', 1, 2)

  dimred1 <- ifelse(dimred[1],1,0)
  dimred2 <- ifelse(dimred[2],1,0)
  dimred3 <- ifelse(dimred[3],1,0)

  xrankfix <- ifelse(is.null(xrankfix),0,xrankfix)
  yrankfix <- ifelse(is.null(yrankfix),0,yrankfix)

  norder <- norder - 1

  fit <- nrrr_cv_my(Y, X, norder, nfold, xr, rfit, xrankfix, yrankfix,
                    jx, jy, p, d, n, maxiter, method,
                    dimred1, dimred2, dimred3, conv, lambda)
  if( sum(fit$rxErrmat)>0 | sum(fit$ryErrmat)>0 | sum(fit$rErrmat)>0 ) stop('Error occurs or the algorithm reaches the maximum iteration')


  return(list(Ag=fit$Ag,Bg=fit$Bg,Al=fit$Al,Bl=fit$Bl,C=fit$C,df=fit$df,
              sse=fit$sse,ic=fit$ic,obj=fit$obj,rx_path=fit$rx_path,ry_path=fit$ry_path,
              r_path=fit$r_path,
              rank=fit$rank,rx=fit$rx,ry=fit$ry,
              rxErrmat=fit$rxErrmat,ryErrmat=fit$ryErrmat,rErrmat=fit$rErrmat))
}
