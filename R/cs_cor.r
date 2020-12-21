# Covariance structure: Compound symmetry
#
# This function generates a covariance matrix with compound symmetry structure.
#
# @param p Dimension.
# @param rho Correlation strength.
# @return A covariance matrix (p by p).
# @export
# @examples
# library(NRRR)
# CorrCS(10, 0.5)
CorrCS <- function(p,rho){
  Sigma <- matrix(nrow=p,ncol=p,rho)
  diag(Sigma) <- 1
  Sigma
}
