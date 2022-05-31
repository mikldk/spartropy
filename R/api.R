strict_convert_integer_matrix <- function(x) {
  stopifnot(is.matrix(x) | is.data.frame(x) | is.vector(x))
  
  x2 <- x
  
  if (!is.data.frame(x2) && !is.matrix(x2)) {
    dim(x2) <- c(length(x2), 1L)
  }
  
  if (!is.integer(x2)) {
    x2 <- df2intmat(x2)
  }
  
  return(x2)
}


#' Entropy using `log()` (natural logarithm)
#' 
#' Se Introduction vignette for definition.
#' 
#' @param d Vector or matrix of observations
#' 
#' @export
entropyE <- function(d) {
  d <- strict_convert_integer_matrix(d)
  ps <- normalise_(frequencies_(d))
  return(entropyE_(ps))
}

#' Entropy using `log2()` 
#' @inherit entropyE
#' @export
entropy2 <- function(d) {
  d <- strict_convert_integer_matrix(d)
  ps <- normalise_(frequencies_(d))
  return(entropy2_(ps))
}

#' Entropy using `log10()`
#' @inherit entropyE
#' @export
entropy10 <- function(d) {
  d <- strict_convert_integer_matrix(d)
  ps <- normalise_(frequencies_(d))
  return(entropy10_(ps))
}

###############################################################

strict_convert_integer_vector <- function(x, max_n) {
  stopifnot(is.vector(x))
  
  x2 <- as.integer(x)
  
  stopifnot(all(abs(x - x2) < 1e-6))
  
  stopifnot(isTRUE(all(x2 >= 1L)))
  stopifnot(isTRUE(all(x2 <= max_n)))
  
  return(x2)
}

strict_convert_integer_Cppindex_vector <- function(x, max_n) {
  x2 <- strict_convert_integer_vector(x, max_n)
  
  # -1 to get 0-based (for C++ indexing)
  x2 <- x2 - 1L
  
  return(x2)
}

strict_convert_integer_Rindex_vector <- function(x, max_n) {
  x2 <- strict_convert_integer_vector(x, max_n)
  
  return(x2)
}

#' Mutual information using `log()` (natural logarithm)
#' 
#' @param d Integer matrix of data where each row is a single observation
#' @param idx_x Indicies (integer vector) for variable set $X$ (in 1 to `ncol(d)`)
#' @param idx_y Indicies (integer vector) for variable set $Y$ (in 1 to `ncol(d)`)
#' 
#' @export
mutinfE <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Cppindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Cppindex_vector(idx_y, ncol(d))
  return(mutual_informationE_implicit_(d, idx_x, idx_y))
}

#' Mutual information using `log2()`
#' @inherit mutinfE
#' @export
mutinf2 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Cppindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Cppindex_vector(idx_y, ncol(d))
  return(mutual_information2_implicit_(d, idx_x, idx_y))
}


#' Mutual information using `log10()`
#' @inherit mutinfE
#' @export
mutinf10 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Cppindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Cppindex_vector(idx_y, ncol(d))
  return(mutual_information10_implicit_(d, idx_x, idx_y))
}

###############################################################

#' Conditional entropy of $X$ given $Y$ using `log()` (natural logarithm)
#' 
#' @inheritParams mutinfE
#' 
#' @export
entropy_condE <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  Hxy <- entropyE(d[, c(idx_x, idx_y)])
  Hy <- entropyE(d[, idx_y])
  
  return(Hxy - Hy)
}

#' Conditional entropy of $X$ given $Y$ using `log2()`
#' @inheritParams mutinfE
#' @export
entropy_cond2 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  Hxy <- entropy2(d[, c(idx_x, idx_y)])
  Hy <- entropy2(d[, idx_y])
  
  return(Hxy - Hy)
}

#' Conditional entropy of $X$ given $Y$ using `log10()`
#' @inheritParams mutinfE
#' @export
entropy_cond10 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  Hxy <- entropy10(d[, c(idx_x, idx_y)])
  Hy <- entropy10(d[, idx_y])
  
  return(Hxy - Hy)
}



###############################################################
infdist_worker <- function(mi, h_x, h_y, h_xy) {
  h_x_y <- h_x - mi
  h_y_x <- h_y - mi
  infdist <- (h_x_y + h_y_x) / h_xy
  return(infdist)
}

#' Shared information distance of $X$ and $Y$ using `log()` (natural logarithm)
#' @inheritParams mutinfE
#' @export
infdistE <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  mi <- mutinfE(d, idx_x, idx_y)
  h_x <- entropyE(d[, idx_x])
  h_y <- entropyE(d[, idx_y])
  h_xy <- entropyE(d[, c(idx_x, idx_y)])

  return(infdist_worker(mi, h_x, h_y, h_xy))
}

#' Shared information distance of $X$ and $Y$ using `log2()` 
#' @inheritParams mutinfE
#' @export
infdist2 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  mi <- mutinf2(d, idx_x, idx_y)
  h_x <- entropy2(d[, idx_x])
  h_y <- entropy2(d[, idx_y])
  h_xy <- entropy2(d[, c(idx_x, idx_y)])
  
  return(infdist_worker(mi, h_x, h_y, h_xy))
}

#' Shared information distance of $X$ and $Y$ using `log10()` 
#' @inheritParams mutinfE
#' @export
infdist10 <- function(d, idx_x, idx_y) {
  d <- strict_convert_integer_matrix(d)
  idx_x <- strict_convert_integer_Rindex_vector(idx_x, ncol(d))
  idx_y <- strict_convert_integer_Rindex_vector(idx_y, ncol(d))
  
  mi <- mutinf10(d, idx_x, idx_y)
  h_x <- entropy10(d[, idx_x])
  h_y <- entropy10(d[, idx_y])
  h_xy <- entropy10(d[, c(idx_x, idx_y)])
  
  return(infdist_worker(mi, h_x, h_y, h_xy))
}
