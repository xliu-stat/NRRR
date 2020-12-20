# Compute || Y - XC ||_F^2
#
# This function computes the sum of squared errors
# with given (A,B,U,V) and (X,Y).
#
# @param Y Response matrix with n rows and jy*d columns.
# @param X Design matrix with n rows and jx*p columns.
# @param Ag Matrix U.
# @param Bg Matrix V.
# @param Al Matrix A.
# @param Bl Matrix B.
# @param jx Number of basis functions to expand x(s).
# @param jy Number of basis functions to expand y(t).
#' @export
# @return The returned items
#   \item{C}{the NRRR estimator of the coefficient matrix C.}
#   \item{XC}{a matrix, prediction of Y.}
#   \item{E}{a matrix, estimated error matrix.}
#   \item{sse,obj}{ a scaler, the sum of squared errors.}
Obj <- function(Y,X,Ag,Bg,Al,Bl,jx,jy){

  if(nrow(Bg)!=ncol(Bg)){
    Bl <- kronecker(diag(jx),Bg)%*%Bl
  }
  if(nrow(Ag)!=ncol(Ag)){
    Al <- kronecker(diag(jy),Ag)%*%Al
  }
  C <- Bl%*%t(Al)
  XC <- X%*%C
  E <- Y - XC
  sse <- sum(E^2)

  return(list(C=C,XC=XC,E=E,sse=sse,obj=sse))
}

