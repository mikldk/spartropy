test_that("frequencies", {
  x <- matrix(1:9, 3, 3)
  
  
  expect_equal(frequencies(x), c(1, 1, 1))
  expect_equal(normalise(frequencies(x)), c(1, 1, 1)/3)

  expect_equal(normalise(c(1, 1, 1)), c(1, 1, 1)/3)
  
  expect_equal(normalise(c(8, 1, 1)), c(8, 1, 1)/10)
  
  ###################
  
  x <- rbind(
    cbind(1, 2, 3),
    cbind(4, 5, 6), 
    cbind(1, 2, 3),
    cbind(4, 5, 6),
    cbind(1, 4, 4),
    cbind(1, 2, 3),
    c(-10, -10, -10)
  )
  storage.mode(x) <- "integer"
  
  expect_equal(sort(frequencies(x)), c(1, 1, 2, 3))
  expect_equal(normalise(sort(frequencies(x))), c(1, 1, 2, 3)/nrow(x))
})

test_that("frequencies_2d", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  storage.mode(x) <- "integer"
  
  expect_equal(frequencies_2d(x, 1, 2), 
               structure(c(3L, 3L, 2L, 2L, 3L, 2L, 3L, 2L, 2L, 1L, 1L, 1L), dim = 4:3, dimnames = list(
                 NULL, c("x", "y", "xy"))))
  
  expect_equal(normalise_2d(frequencies_2d(x, 1, 2)), 
               structure(c(0.6, 0.6, 0.4, 0.4, 0.6, 0.4, 0.6, 0.4, 0.4, 0.2, 
                           0.2, 0.2), dim = 4:3))
})

