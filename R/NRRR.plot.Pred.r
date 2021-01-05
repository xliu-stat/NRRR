#' @title
#' Plot the predicted response trajectory
#'
#' @description
#' This function plots the predicted response trajectory based on the output of \code{\link{NRRR.pred}}.
#'
#'
#' @usage
#' NRRR.plot.pred(Ypred, Y = NULL, i_ind = 1, yi_ind = 1,
#'                tseq, t_index = NULL, x_lab = NULL,
#'                y_lab = NULL)
#'
#'
#' @param Ypred an array of dimension \code{(n, d, length(tseq))} where n is
#'              the sample size and d is the number of components in the
#'              multivariate response. The response trajectory is predicted at
#'              a set of time points \code{tseq}.
#' @param Y an array of dimension \code{(n, d, length(tseq))}. It is the truly
#'          observed response trajectory in the form of discrete observations at
#'          time points \code{tseq}. This term is only available
#'          when the function is applied to a testing set.
#' @param i_ind a user-specified sample ID for which the plot is to be drawn.
#'              If the prediction is only conducted on one sample,
#'              then set \code{i_ind = 1}; otherwise it should be an integer
#'              less than or equal to n. Default is 1.
#' @param yi_ind a user-specified response component index for which the plot
#'               is to be drawn. It should satisfy \eqn{0 < yi_ind \le d}.
#'               Default is 1.
#' @param tseq a sequence of time points at which the predicted response values
#'             are obtained.
#' @param x_lab,y_lab the user-specified x-axis and y-axis label,
#'                    and it should be given as a character string, e.g., x_lab = "Time".
#'                    Default is NULL.
#' @param t_index the user-specified x-axis tick marks, and it should be
#'                given as a vector of character strings of the same length as tseq.
#'                Default is NULL.
#'
#' @return A ggplot2 object. If the truly observed response trajectory is available, then
#'         it is plotted as a black line while the predicted trajectory is in red.
#'
#'
#' @details
#' An example of its usage can be found in the vignette of electricity demand analysis.
#'
#' @references
#' Liu, X., Ma, S., & Chen, K. (2020). Multivariate Functional Regression via Nested Reduced-Rank Regularization.
#' arXiv: Methodology.
#'
#' @import ggplot2
#' @export

NRRR.plot.pred <- function(Ypred, Y = NULL, i_ind = 1, yi_ind = 1,
                           tseq, t_index = NULL, x_lab = NULL,
                           y_lab = NULL){

  if (i_ind > dim(Ypred)[1]) stop("i_ind cannot be greater than n")
  if (yi_ind > dim(Ypred)[2]) stop("yi_ind cannot be greater than d")

  iy <- i_ind
  l <- yi_ind

  tseq <- factor(tseq)
  if (!is.null(t_index)) {
    levels(tseq) <- t_index
  } else {
    t_index <- tseq
  }

  if (is.null(x_lab)) x_lab <- "t"
  if (is.null(y_lab)) y_lab <- "Predicted Value"

  if (is.null(Y)) {
    b <- Ypred[iy,l,]
    y.range <- range(b)
    cur <- data.frame(b,tseq)
    ggplot(data = cur, aes(x = tseq)) +
      geom_line(aes(y = b,group = 1), colour = "red", linetype = 2,size=1.5) +
      scale_x_discrete(breaks = t_index[seq(1,length(t_index),2)]) +
      xlab(x_lab) +
      ylab(y_lab) +
      ylim(y.range) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.text=element_text(size=12),axis.title=element_text(size=20),
            plot.title = element_text(hjust = 0.5, size = 20))
  } else {
    a <- Y[iy,l,]
    b <- Ypred[iy,l,]
    y.range <- range(a,b)
    cur <- data.frame(a,b,tseq)
    ggplot(data = cur, aes(x = tseq)) +
      geom_line(aes(y = a,group = 1), colour = "black",size=1.5) +
      geom_line(aes(y = b,group = 1), colour = "red", linetype = 2,size=1.5) +
      scale_x_discrete(breaks = t_index[seq(1,length(t_index),2)]) +
      xlab(x_lab) +
      ylab(y_lab) +
      ylim(y.range) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(axis.text=element_text(size=12),axis.title=element_text(size=20),
            plot.title = element_text(hjust = 0.5, size = 20))
  }
}


