test_that("entropy", {
  expect_equal(entropy(c(0.5, 0.5)), 0.693147180559945)
  expect_equal(entropy2(c(0.5, 0.5)), 1)
  expect_equal(entropy10(c(0.5, 0.5)), 0.301029995663981)

  x <- 1:10
  p <- x / sum(x)
  expect_equal(entropy(p), 2.15128172065184)
  
  x <- c(NA, 1:10)
  p <- x / sum(x)
  expect_true(is.na(entropy(p)))
  
  
  # http://www.cs.tau.ac.il/~iftachh/Courses/Info/Fall14/Printouts/Lesson2_h.pdf
  expect_equal(entropy2(c(1/2, 1/2, 1/4)), 1.5)
})




test_that("mutual_information", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  
  ps <- frequencies_2d(x, 0, 1) |> normalise_2d()
  mi <- mutual_information2(ps)
  
  
  H_x <- entropy2(normalise(frequencies(x[, 1, drop = FALSE])))
  H_y <- entropy2(normalise(frequencies(x[, 2, drop = FALSE])))
  H_xy <- entropy2(normalise(frequencies(x)))
  
  H_y_x <- H_xy - H_x
  H_x_y <- H_xy - H_y
  
  expect_equal(H_y_x, H_x_y)
  
  # mi
  # H_x - H_x_y
  # H_y - H_y_x
  # H_x + H_y - H_xy
  # H_xy - H_x_y - H_y_x
  expect_equal(mi, H_x - H_x_y)
  expect_equal(mi, H_y - H_y_x)
  expect_equal(mi, H_x + H_y - H_xy)
  expect_equal(mi, H_xy - H_x_y - H_y_x)
})


test_that("mutual_information", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  
  ps <- frequencies_2d(x, 0, 1) |> normalise_2d()
  mi <- mutual_information2(ps)
  
  expect_equal(mi, mutual_information2_implicit(x, 0, 1))
  
  ####
  
  x2 <- rbind(x, x, x, x, x, x, x)
  x3 <- rbind(x2, x2 + 1, x2 + 1, x2 + 1, x2 + 2)
  x4 <- rbind(x3, x3, x3 + 4, x3 + 8, x3 + 10, x3, x3)
  x5 <- rbind(x4, x4, x4, x4, x4, x4, x4)
  x6 <- cbind(x5, rev(x5[, 1]))
  x_big <- x6
  
  expect_equal(frequencies_2d(x_big, 0, 1) |> normalise_2d() |> mutual_information2(), 
               mutual_information2_implicit(x_big, 0, 1))
  
  expect_equal(frequencies_2d(x_big, c(0, 2), 1) |> normalise_2d() |> mutual_information2(), 
               mutual_information2_implicit(x_big, c(0, 2), 1))
  
  expect_equal(frequencies_2d(x_big, c(0, 1), 2) |> normalise_2d() |> mutual_information2(), 
               mutual_information2_implicit(x_big, c(0, 1), 2))
  
  expect_equal(frequencies_2d(x_big, 0, c(1, 2)) |> normalise_2d() |> mutual_information2(), 
               mutual_information2_implicit(x_big, 0, c(1, 2)))
  
  if (FALSE) {
    microbenchmark::microbenchmark(
      frequencies_2d(x_big, c(0, 2), 1) |> normalise_2d() |> mutual_information2(),
      mutual_information2_implicit(x_big, c(0, 2), 1)
    )
  }
})


