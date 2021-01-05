#' @title
#' Plot heatmap for the functional regression surface
#'
#' @description
#' This function creates heatmaps for the functional regression surface in a
#' multivariate functional linear regression. Based on the fitting results from the
#' nested reduced-rank regression, different kinds of regression surfaces
#' (at the original scale or the latent scale) can be visualized to give a
#' clear illustration of the functional correlation between the user-specified
#' predictor (or latent predictor) trajectory and response
#' (or latent response) trajectory.
#'
#'
#' @usage
#' NRRR.plot.reg(Ag, Bg, Al, Bl, rx, ry, sseq, phi, tseq, psi,
#'               x_ind, y_ind, x_lab = NULL, y_lab = NULL,
#'               tseq_index = NULL, sseq_index = NULL,
#'               method = c("latent", "x_original",
#'               "y_original", "original")[1])
#'
#'
#' @param Ag,Bg,Al,Bl,rx,ry the estimated U, V, A, B, rx and ry from a NRRR fitting.
#' @param sseq the sequence of time points at which the predictor trajectory is observed.
#' @param phi the set of basis functions to expand the predictor trajectory.
#' @param tseq the sequence of time points at which the response trajectory is observed.
#' @param psi the set of basis functions to expand the response trajectory.
#' @param x_ind,y_ind two indices to locate the regression surface for which the heat map is to be drawn.
#'                    If \code{method = "original"}, then \eqn{0 < x_ind <= p, 0 < y_ind <= d}
#'                    and the function plots \eqn{C_{x_ind,y_ind}(s,t)} in Eq. (1) of the NRRR paper.
#'                    If \code{method = "latent"}, then \eqn{0 < x_ind <= rx, 0 < y_ind <= ry}
#'                    and the function plots \eqn{C^*_{x_ind,y_ind}(s,t)} in Eq. (2) of the NRRR paper.
#'                    If \code{method = "y_original"}, then \eqn{0 < x_ind <= rx, 0 < y_ind <= d}.
#'                    If \code{method = "x_original"}, then \eqn{0 < x_ind <= p, 0 < y_ind <= ry}.
#' @param x_lab,y_lab the user-specified x-axis (with x_lab for predictor) and
#'                    y-axis (with y_lab for response) label,
#'                    and it should be given as a character string, e.g., x_lab = "Temperature".
#' @param tseq_index,sseq_index the user-specified x-axis (with sseq_index for predictor)
#'                              and y-axis (with tseq_index for response) tick marks, and it should be
#'                              given as a vector of character strings of the same length as sseq or tseq, respectively.
#' @param method 'original': the function plots the correlation heatmap between the original
#'               functional response \eqn{y_i(t)} and the original functional predictor \eqn{x_j(s)};
#'               'latent': the function plots the correlation heatmap between
#'               the latent functional response \eqn{y^*_i(t)} and the latent functional predictor \eqn{x^*_j(s)};
#'               'y_original': the function plots the correlation heatmap between \eqn{y_i(t)} and \eqn{x^*_j(s)};
#'               'x_original': the function plots the correlation heatmap between \eqn{y^*_i(t)} and \eqn{x_j(s)}.
#'
#' @return A ggplot2 object.
#'
#' @details
#' More details and the examples of its usage can be found in the vignette of electricity demand analysis.
#'
#' @references
#' Liu, X., Ma, S., & Chen, K. (2020). Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export

NRRR.plot.reg <- function(Ag, Bg, Al, Bl, rx, ry,
                                 sseq, phi, tseq, psi,
                                 x_ind, y_ind,
                                 x_lab = NULL, y_lab = NULL,
                                 tseq_index = NULL, sseq_index = NULL,
                                 method = c("latent", "x_original", "y_original","original")[1]
){
  ry <- ry
  rx <- rx
  Al <- Al
  Bl <- Bl
  Ag <- Ag
  Bg <- Bg
  p <- dim(Bg)[1]
  d <- dim(Ag)[1]
  ns <- dim(phi)[1]
  nt <- dim(psi)[1]
  jx <- dim(phi)[2]
  jy <- dim(psi)[2]

  if (method == "original" & any(c(0 > x_ind, x_ind > p, 0 > y_ind, y_ind > d))) stop("when 'original' is selected, 0 < x_ind <= p, 0 < y_ind <= d")
  if (method == "latent" & any(c(0 > x_ind, x_ind > rx, 0 > y_ind, y_ind > ry))) stop("when 'latent' is selected, 0 < x_ind <= rx, 0 < y_ind <= ry")
  if (method == "y_original" & any(c(0 > x_ind, x_ind > rx, 0 > y_ind, y_ind > d))) stop("when 'y_original' is selected, 0 < x_ind <= rx, 0 < y_ind <= d")
  if (method == "x_original" & any(c(0 > x_ind, x_ind > p, 0 > y_ind, y_ind > ry))) stop("when 'x_original' is selected, 0 < x_ind <= p, 0 < y_ind <= ry")


  Jpsi <- matrix(nrow = jy, ncol = jy, 0)
  tdiff <- (tseq - c(0, tseq[-nt]))
  for (t in 1:nt) Jpsi <- Jpsi + psi[t, ] %*% t(psi[t, ]) * tdiff[t]
  eJpsi <- eigen(Jpsi)
  Jpsihalf <- eJpsi$vectors %*% diag(sqrt(eJpsi$values)) %*% t(eJpsi$vectors)
  Jpsihalfinv <- eJpsi$vectors %*% diag(1 / sqrt(eJpsi$values)) %*% t(eJpsi$vectors)


  alindex <- rep(1:ry,jy)
  Alstar <- Al[order(alindex),]
  blindex <- rep(1:rx,jx)
  Blstar <- Bl[order(blindex),]


  Alstar <- kronecker(diag(ry),Jpsihalfinv)%*%Alstar
  Core <- Alstar%*%t(Blstar)


  comp <- array(dim = c(ry, rx, nt, ns), NA)

  for (i in 1:ry){
    for (j in 1:rx){
      comp[i,j,,] <- psi%*%Core[c(((i-1)*jy + 1):(i*jy)),
                                    c(((j-1)*jx + 1):(j*jx))]%*%t(phi) # nt \times ns
    }
  }


  if (method == "original"){
    comp.ori <- array(dim = c(d, p, nt, ns), NA)
    for (i in 1:nt){
      for (j in 1:ns){
        comp.ori[ , , i, j] <- Ag%*%comp[ , , i, j]%*%t(Bg)
      }
    }
    comp <- comp.ori
  } else if (method == "y_original"){
    comp.ori <- array(dim = c(d, rx, nt, ns), NA)
    for (i in 1:nt){
      for (j in 1:ns){
        comp.ori[ , , i, j] <- Ag%*%comp[ , , i, j]
      }
    }
    comp <- comp.ori
  } else if (method == "x_original"){
    comp.ori <- array(dim = c(ry, p, nt, ns), NA)
    for (i in 1:nt){
      for (j in 1:ns){
        comp.ori[ , , i, j] <- comp[ , , i, j]%*%t(Bg)
      }
    }
    comp <- comp.ori
  }

  # heatmap for coefficient function
  # library(ggplot2)
  # library(reshape2)

    dat <- reshape2::melt(as.matrix(comp[y_ind, x_ind, , ]))

    dat$Var1 <- factor(dat$Var1)
    dat$Var2 <- factor(dat$Var2)
    if (is.null(tseq_index)) {
      levels(dat$Var1) <- factor(tseq)
      y_breaks <- tseq[seq(1,nt,2)]
    } else {
      levels(dat$Var1) <- tseq_index
      y_breaks <- tseq_index[seq(1,nt,2)]
    }

    if (is.null(sseq_index)) {
      levels(dat$Var2) <- factor(sseq)
      x_breaks <- sseq[seq(1,ns,2)]
    } else {
      levels(dat$Var2) <- sseq_index
      x_breaks <- sseq_index[seq(1,ns,2)]
    }

    # Var1 <- dat$Var1
    # Var2 <- dat$Var2
    # value <- dat$value


    ggplot2::ggplot(dat, aes(Var2, Var1, fill = value)) + geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white",
                           high = "red", midpoint = 0, limits= range(comp),
                           space = "Lab", name = "",
                           na.value = "grey50", guide = "colourbar", aesthetics = "fill")+
      scale_x_discrete(breaks = x_breaks) +
      scale_y_discrete(breaks = y_breaks) +
      ylab(ifelse(is.null(y_lab), "response (t)", y_lab)) +
      xlab(ifelse(is.null(x_lab), "predictor (s)", x_lab)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.text=element_text(size=10),axis.title=element_text(size=13),
            plot.title = element_text(hjust = 0.5, size = 15))
}
