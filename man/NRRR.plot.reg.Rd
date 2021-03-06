% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NRRR.plot.RegSurface.r
\name{NRRR.plot.reg}
\alias{NRRR.plot.reg}
\title{Plot heatmap for the functional regression surface}
\usage{
NRRR.plot.reg(Ag, Bg, Al, Bl, rx, ry, sseq, phi, tseq, psi,
              x_ind, y_ind, x_lab = NULL, y_lab = NULL,
              tseq_index = NULL, sseq_index = NULL,
              method = c("latent", "x_original",
              "y_original", "original")[1])
}
\arguments{
\item{Ag, Bg, Al, Bl, rx, ry}{the estimated U, V, A, B, rx and ry from a NRRR fitting.}

\item{sseq}{the sequence of time points at which the predictor trajectory is observed.}

\item{phi}{the set of basis functions to expand the predictor trajectory.}

\item{tseq}{the sequence of time points at which the response trajectory is observed.}

\item{psi}{the set of basis functions to expand the response trajectory.}

\item{x_ind, y_ind}{two indices to locate the regression surface for which the heat map is to be drawn.
If \code{method = "original"}, then \eqn{0 < x_ind <= p, 0 < y_ind <= d}
and the function plots \eqn{C_{x_ind,y_ind}(s,t)} in Eq. (1) of the NRRR paper.
If \code{method = "latent"}, then \eqn{0 < x_ind <= rx, 0 < y_ind <= ry}
and the function plots \eqn{C^*_{x_ind,y_ind}(s,t)} in Eq. (2) of the NRRR paper.
If \code{method = "y_original"}, then \eqn{0 < x_ind <= rx, 0 < y_ind <= d}.
If \code{method = "x_original"}, then \eqn{0 < x_ind <= p, 0 < y_ind <= ry}.}

\item{x_lab, y_lab}{the user-specified x-axis (with x_lab for predictor) and
y-axis (with y_lab for response) label,
and it should be given as a character string, e.g., x_lab = "Temperature".}

\item{tseq_index, sseq_index}{the user-specified x-axis (with sseq_index for predictor)
and y-axis (with tseq_index for response) tick marks, and it should be
given as a vector of character strings of the same length as sseq or tseq, respectively.}

\item{method}{'original': the function plots the correlation heatmap between the original
functional response \eqn{y_i(t)} and the original functional predictor \eqn{x_j(s)};
'latent': the function plots the correlation heatmap between
the latent functional response \eqn{y^*_i(t)} and the latent functional predictor \eqn{x^*_j(s)};
'y_original': the function plots the correlation heatmap between \eqn{y_i(t)} and \eqn{x^*_j(s)};
'x_original': the function plots the correlation heatmap between \eqn{y^*_i(t)} and \eqn{x_j(s)}.}
}
\value{
A ggplot2 object.
}
\description{
This function creates heatmaps for the functional regression surface in a
multivariate functional linear regression. Based on the fitting results from the
nested reduced-rank regression, different kinds of regression surfaces
(at the original scale or the latent scale) can be visualized to give a
clear illustration of the functional correlation between the user-specified
predictor (or latent predictor) trajectory and response
(or latent response) trajectory.
}
\details{
More details and the examples of its usage can be found in the vignette of electricity demand analysis.
}
\references{
Liu, X., Ma, S., & Chen, K. (2020). Multivariate Functional Regression via Nested Reduced-Rank Regularization.
arXiv: Methodology.
}
