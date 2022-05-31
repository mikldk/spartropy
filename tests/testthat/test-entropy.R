test_that("entropy", {
  expect_equal(entropyE_(c(0.5, 0.5)), 0.693147180559945)
  expect_equal(entropy2_(c(0.5, 0.5)), 1)
  expect_equal(entropy10_(c(0.5, 0.5)), 0.301029995663981)

  x <- 1:10
  p <- x / sum(x)
  expect_equal(entropyE_(p), 2.15128172065184)
  
  x <- c(NA, 1:10)
  p <- x / sum(x)
  expect_true(is.na(entropyE_(p)))
  
  
  # http://www.cs.tau.ac.il/~iftachh/Courses/Info/Fall14/Printouts/Lesson2_h.pdf
  expect_equal(entropy2_(c(1/2, 1/2, 1/4)), 1.5)
})




test_that("mutual_information", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  
  ps <- normalise_2d_(frequencies_2d_(x, 0, 1))
  mi1 <- mutual_information2_(ps)
  mi2 <- mutual_information2_implicit_(x, 0, 1)
  expect_equal(mi1, mi2)
  
  H_x <- entropy2_(normalise_(frequencies_(x[, 1, drop = FALSE])))
  H_y <- entropy2_(normalise_(frequencies_(x[, 2, drop = FALSE])))
  H_xy <- entropy2_(normalise_(frequencies_(x)))
  
  H_y_x <- H_xy - H_x
  H_x_y <- H_xy - H_y
  
  expect_equal(H_y_x, H_x_y)
  
  # mi
  # H_x - H_x_y
  # H_y - H_y_x
  # H_x + H_y - H_xy
  # H_xy - H_x_y - H_y_x
  expect_equal(mi1, H_x - H_x_y)
  expect_equal(mi1, H_y - H_y_x)
  expect_equal(mi1, H_x + H_y - H_xy)
  expect_equal(mi1, H_xy - H_x_y - H_y_x)
  
  
  ############################
  
  miEa <- mutual_informationE_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  miEb <- mutual_informationE_implicit_(x, 0, 1)
  expect_equal(miEa, miEb)
  expect_equal(miEa, miEb)
  
  mi2a <- mutual_information2_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  mi2b <- mutual_information2_implicit_(x, 0, 1)
  expect_equal(mi2a, mi2b)
  
  mi10a <- mutual_information10_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  mi10b <- mutual_information10_implicit_(x, 0, 1)
  expect_equal(mi10a, mi10b)
})


test_that("mutual_information", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  
  ps <- normalise_2d_(frequencies_2d_(x, 0, 1))
  mi <- mutual_information2_(ps)
  
  expect_equal(mi, mutual_information2_implicit_(x, 0, 1))
  
  ####
  
  x2 <- rbind(x, x, x, x, x, x, x)
  x3 <- rbind(x2, x2 + 1, x2 + 1, x2 + 1, x2 + 2)
  x4 <- rbind(x3, x3, x3 + 4, x3 + 8, x3 + 10, x3, x3)
  x5 <- rbind(x4, x4, x4, x4, x4, x4, x4)
  x6 <- cbind(x5, rev(x5[, 1]))
  x_big <- x6
  storage.mode(x_big) <- "integer"
  
  
  # E:
  expect_equal(mutual_informationE_(normalise_2d_(frequencies_2d_(x_big, 0, 1))), 
               mutual_informationE_implicit_(x_big, 0, 1))
  
  expect_equal(mutual_informationE_(normalise_2d_(frequencies_2d_(x_big, c(0, 2), 1))),
               mutual_informationE_implicit_(x_big, c(0, 2), 1))
  
  expect_equal(mutual_informationE_(normalise_2d_(frequencies_2d_(x_big, c(0, 1), 2))),
               mutual_informationE_implicit_(x_big, c(0, 1), 2))
  
  expect_equal(mutual_informationE_(normalise_2d_(frequencies_2d_(x_big, 0, c(1, 2)))),
               mutual_informationE_implicit_(x_big, 0, c(1, 2)))
  
  # 2:
  expect_equal(mutual_information2_(normalise_2d_(frequencies_2d_(x_big, 0, 1))), 
               mutual_information2_implicit_(x_big, 0, 1))
  
  expect_equal(mutual_information2_(normalise_2d_(frequencies_2d_(x_big, c(0, 2), 1))),
               mutual_information2_implicit_(x_big, c(0, 2), 1))
  
  expect_equal(mutual_information2_(normalise_2d_(frequencies_2d_(x_big, c(0, 1), 2))),
               mutual_information2_implicit_(x_big, c(0, 1), 2))
  
  expect_equal(mutual_information2_(normalise_2d_(frequencies_2d_(x_big, 0, c(1, 2)))),
               mutual_information2_implicit_(x_big, 0, c(1, 2)))
  
  # 10:
  expect_equal(mutual_information10_(normalise_2d_(frequencies_2d_(x_big, 0, 1))), 
               mutual_information10_implicit_(x_big, 0, 1))
  
  expect_equal(mutual_information10_(normalise_2d_(frequencies_2d_(x_big, c(0, 2), 1))),
               mutual_information10_implicit_(x_big, c(0, 2), 1))
  
  expect_equal(mutual_information10_(normalise_2d_(frequencies_2d_(x_big, c(0, 1), 2))),
               mutual_information10_implicit_(x_big, c(0, 1), 2))
  
  expect_equal(mutual_information10_(normalise_2d_(frequencies_2d_(x_big, 0, c(1, 2)))),
               mutual_information10_implicit_(x_big, 0, c(1, 2)))
  
  #
  if (FALSE) {
    microbenchmark::microbenchmark(
      mutual_information2_(normalise_2d_(frequencies_2d(x_big, c(1, 3), 2))),
      mutual_information2_implicit_(x_big, c(1, 3), 2)
    )
  }
})



test_that("entropy", {
  x <- df2intmat(mtcars)
  n <- frequencies_(x)
  e <- entropyE_(normalise_(n))
  
  expect_equal(e, 3.46573590279973)
})

test_that("entropy package", {
  if (!require("entropy", quietly = TRUE)) {
    skip(message = "Did not have 'entropy' package -- skipping comparisons")
  }
  
  x <- df2intmat(mtcars[, 1:2])
  n <- frequencies_(x)
  
  e1 <- entropy::entropy.empirical(n, unit = "log2")
  e2 <- entropy2_(normalise_(n))
  
  expect_equal(e1, e2)
  
  ##
  
  y2d <- table(x[, 1], x[, 2])
  
  mi1 <- entropy::mi.empirical(y2d, unit = "log2")
  mi2 <- mutual_information2_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  expect_equal(mi1, mi2)
  mi2 <- mutual_information2_implicit_(x, 0, 1)
  expect_equal(mi1, mi2)
  
  mi1 <- entropy::mi.empirical(y2d, unit = "log")
  mi2 <- mutual_informationE_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  expect_equal(mi1, mi2)
  mi2 <- mutual_informationE_implicit_(x, 0, 1)
  expect_equal(mi1, mi2)
  
  mi1 <- entropy::mi.empirical(y2d, unit = "log10")
  mi2 <- mutual_information10_(normalise_2d_(frequencies_2d_(x, 0, 1)))
  expect_equal(mi1, mi2)
  mi2 <- mutual_information10_implicit_(x, 0, 1)
  expect_equal(mi1, mi2)
})


test_that("infotheo package", {
  if (!require("infotheo", quietly = TRUE)) {
    skip(message = "Did not have 'infotheo' package -- skipping comparisons")
  }
  
  x <- df2intmat(mtcars[, c(1, 4, 6)])
  
  #####
  
  mi1 <- infotheo::mutinformation(X = x[, 1], Y = x[, c(2, 3)], method = "emp")
  mi2 <- mutinfE(mtcars, 1, c(4, 6))
  expect_equal(mi1, mi2)
  
  mi1 <- infotheo::mutinformation(X = x[, c(1, 2)], Y = x[, 3], method = "emp")
  mi2 <- mutinfE(mtcars, c(1, 4), 6)
  expect_equal(mi1, mi2)
  
  mi1 <- infotheo::mutinformation(X = x[, c(1, 3)], Y = x[, 2], method = "emp")
  mi2 <- mutinfE(mtcars, c(1, 6), 4)
  expect_equal(mi1, mi2)
  
  #####
  
  ce1 <- infotheo::condentropy(X = x[, 2], Y = x[, 1], method = "emp")
  ce2 <- entropy_condE(mtcars, 4, 1)
  expect_equal(ce1, ce2)
})

