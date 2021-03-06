% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api.R
\name{mutinf10}
\alias{mutinf10}
\title{Mutual information using `log10()`}
\usage{
mutinf10(d, idx_x, idx_y)
}
\arguments{
\item{d}{Integer matrix of data where each row is a single observation}

\item{idx_x}{Indicies (integer vector) for variable set $X$ (in 1 to `ncol(d)`)}

\item{idx_y}{Indicies (integer vector) for variable set $Y$ (in 1 to `ncol(d)`)}
}
\description{
Mutual information using `log10()`
}
