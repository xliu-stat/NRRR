# Covariance structure: Autoregressive
#
# This function generates a covariance matrix with autoregressive structure.
#
# @param p Dimension.
# @param rho Correlation strength.
# @return A covariance matrix (p by p).
# @export
# @examples
# library(NRRR)
# CorrAR(10, 0.5)
CorrAR <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,NA)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  Sigma
}

