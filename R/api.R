require_integer_matrix <- function(x) {
  stopifnot(is.matrix(x))
  stopifnot(is.integer(x))
}

strict_convert_integer_vector <- function(x) {
  stopifnot(is.vector(x))
  
  x2 <- as.integer(x)
  
  stopifnot(all(abs(x - x2) < 1e-6))
  
  return(x2)
}

##############################

#' Entropy using `log()` (natural logarithm)
#' 
#' @param p Vector of probabilities that constitutes a distribution (all >= 0 and <= 1 and sums to 1)
#' 
#' @export
entropy <- function(p) {
  return(entropy_(p))
}

#' Entropy using `log2()` 
#' @inheritParams entropy
#' @export
entropy2 <- function(p) {
  return(entropy2_(p))
}

#' Entropy using `log10()` 
#' @inheritParams entropy
#' @export
entropy10 <- function(p) {
  return(entropy10_(p))
}

##############################


#' Frequencies
#' 
#' Similar to `table()`
#' 
#' @param x Integer matrix of data where each row is a single observation
#' 
#' @export
frequencies <- function(x) {
  require_integer_matrix(x)
  
  return(frequencies_(x))
}

#' Normalise counts
#' 
#' @param x Integer vector of counts from [frequencies()]
#' 
#' @export
normalise <- function(x) {
  x <- strict_convert_integer_vector(x)
  
  return(normalise_(x))
}


#' Frequencies
#' 
#' Similar to `table()` for two sets of variables, $X$ and $Y$, 
#' but returns a matrix with counts instead
#' 
#' @param x Integer matrix of data where each row is a single observation
#' @param idx_x Indicies (integer vector) for variable set $X$
#' @param idx_y Indicies (integer vector) for variable set $Y$
#' 
#' @export
frequencies_2d <- function(x, idx_x, idx_y) {
  require_integer_matrix(x)
  
  idx_x <- strict_convert_integer_vector(idx_x)
  idx_y <- strict_convert_integer_vector(idx_y)

  stopifnot(isTRUE(all(idx_x >= 1L)))
  stopifnot(isTRUE(all(idx_y >= 1L)))
  stopifnot(isTRUE(all(idx_x <= ncol(x))))
  stopifnot(isTRUE(all(idx_y <= ncol(x))))
  
  # -1 to get 0-based
  y <- frequencies_2d_(x, idx_x - 1L, idx_y - 1L)
  
  return(y)
}

#' Normalise counts
#' 
#' @param x Integer matrix of counts from [frequencies_2d()]
#' 
#' @export
normalise_2d <- function(x) {
  require_integer_matrix(x)
  
  return(normalise_2d_(x))
}

####################

#' Mutual information using `log()` (natural logarithm)
#' 
#' @param Matrix of probabilities obtained from [frequencies_2d()] and [normalise_2d()]
#' 
#' @export
mutual_information <- function(ps) {
  y <- mutual_information_(ps)
  return(y)
}

#' Mutual information using `log2()`
#' 
#' @inheritParams mutual_information
#' 
#' @export
mutual_information2 <- function(ps) {
  y <- mutual_information2_(ps)
  return(y)
}

#' Mutual information using `log10()`
#' 
#' @inheritParams mutual_information
#' 
#' @export
mutual_information10 <- function(ps) {
  y <- mutual_information10_(ps)
  return(y)
}


####################

#' Mutual information using `log()` (natural logarithm)
#' 
#' @inheritParams frequencies_2d
#' 
#' @export
mutual_information_implicit <- function(x, idx_x, idx_y) {
  require_integer_matrix(x)
  
  idx_x <- strict_convert_integer_vector(idx_x)
  idx_y <- strict_convert_integer_vector(idx_y)
  
  stopifnot(isTRUE(all(idx_x >= 1L)))
  stopifnot(isTRUE(all(idx_y >= 1L)))
  stopifnot(isTRUE(all(idx_x <= ncol(x))))
  stopifnot(isTRUE(all(idx_y <= ncol(x))))
  
  # -1 to get 0-based
  y <- mutual_information_implicit_(x, idx_x - 1L, idx_y - 1L)
  
  return(y)
}

#' Mutual information using `log2()`
#' 
#' @inheritParams frequencies_2d
#' 
#' @export
mutual_information2_implicit <- function(x, idx_x, idx_y) {
  require_integer_matrix(x)
  
  idx_x <- strict_convert_integer_vector(idx_x)
  idx_y <- strict_convert_integer_vector(idx_y)
  
  stopifnot(isTRUE(all(idx_x >= 1L)))
  stopifnot(isTRUE(all(idx_y >= 1L)))
  stopifnot(isTRUE(all(idx_x <= ncol(x))))
  stopifnot(isTRUE(all(idx_y <= ncol(x))))
  
  # -1 to get 0-based
  y <- mutual_information2_implicit_(x, idx_x - 1L, idx_y - 1L)
  
  return(y)
}

#' Mutual information using `log10()`
#' 
#' @inheritParams frequencies_2d
#' 
#' @export
mutual_information10_implicit <- function(x, idx_x, idx_y) {
  require_integer_matrix(x)
  
  idx_x <- strict_convert_integer_vector(idx_x)
  idx_y <- strict_convert_integer_vector(idx_y)
  
  stopifnot(isTRUE(all(idx_x >= 1L)))
  stopifnot(isTRUE(all(idx_y >= 1L)))
  stopifnot(isTRUE(all(idx_x <= ncol(x))))
  stopifnot(isTRUE(all(idx_y <= ncol(x))))
  
  # -1 to get 0-based
  y <- mutual_information10_implicit_(x, idx_x - 1L, idx_y - 1L)
  
  return(y)
}
