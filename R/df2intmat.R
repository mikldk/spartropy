#' Converts data frame to integer matrix
#' 
#' The conversion is done by converting to characters, then a factor and then integer
#' 
#' @examples 
#' x <- mtcars[1:4, 1:5]
#' x
#' df2intmat(x)
#' 
#' @param x Data frame
df2intmat <- function(x) {
  y <- apply(x, 2, function(z) as.integer(factor(as.character(z))))
  return(y)
}
