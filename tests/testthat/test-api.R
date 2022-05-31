test_that("strict_convert_integer_matrix", {
  x <- strict_convert_integer_matrix(mtcars[1:4, 1:4])
  expect_true(is.matrix(x))
  expect_true(is.integer(x))
  expect_equal(x, structure(c(1L, 1L, 3L, 2L, 2L, 2L, 1L, 2L, 
                              2L, 2L, 1L, 3L, 1L, 1L, 2L, 1L), 
                            dim = c(4L, 4L), 
                            dimnames = list(NULL, c("mpg", "cyl", "disp", "hp"))))
  
  #################
  
  x <- strict_convert_integer_matrix(mtcars[1:4, 1])
  expect_true(is.matrix(x))
  expect_true(is.integer(x))
  expect_equal(x, structure(c(1L, 1L, 3L, 2L), dim = c(4L, 1L)))
})

test_that("strict_convert_integer_vector", {
  expect_error(strict_convert_integer_vector(1:10, 5))
  expect_error(strict_convert_integer_Cppindex_vector(1:10, 5))
  expect_error(strict_convert_integer_Rindex_vector(1:10, 5))
  
  expect_equal(strict_convert_integer_vector(1:3, 3), 1:3)
  expect_equal(strict_convert_integer_Cppindex_vector(1:3, 3), 0:2)
  expect_equal(strict_convert_integer_Rindex_vector(1:3, 3), 1:3)
})

test_that("API", {
  x <- matrix(c(
    1, 1, 
    2, 1, 
    1, 2, 
    2, 2, 
    2, 2
  ), ncol = 2, byrow = TRUE)
  
  expect_equal(entropyE(x[, 1]), entropyE_(normalise_(frequencies_(x[, 1, drop = FALSE]))))
  expect_equal(entropy2(x[, 1]), entropy2_(normalise_(frequencies_(x[, 1, drop = FALSE]))))
  expect_equal(entropy10(x[, 1]), entropy10_(normalise_(frequencies_(x[, 1, drop = FALSE]))))
  
  
  expect_equal(entropyE(x[, c(1, 2)]), entropyE(x[, c(2, 1)]))
  expect_equal(entropyE(x[, c(1, 2)]), entropyE_(normalise_(frequencies_(x[, c(1, 2), drop = FALSE]))))
  
  expect_equal(entropy2(x[, c(1, 2)]), entropy2(x[, c(2, 1)]))
  expect_equal(entropy2(x[, c(1, 2)]), entropy2_(normalise_(frequencies_(x[, c(1, 2), drop = FALSE]))))
  
  expect_equal(entropy10(x[, c(1, 2)]), entropy10(x[, c(2, 1)]))
  expect_equal(entropy10(x[, c(1, 2)]), entropy10_(normalise_(frequencies_(x[, c(1, 2), drop = FALSE]))))
  
  ###########
  
  # Note R-indexing vs C++-indexing
  expect_equal(mutinfE(x, 1, 2), mutual_informationE_(normalise_2d_(frequencies_2d_(x, 0, 1))))
  expect_equal(mutinfE(x, 1, 2), mutinfE(x, 2, 1))
  
  expect_equal(mutinf2(x, 1, 2), mutual_information2_(normalise_2d_(frequencies_2d_(x, 0, 1))))
  expect_equal(mutinf2(x, 1, 2), mutinf2(x, 2, 1))
  
  expect_equal(mutinf10(x, 1, 2), mutual_information10_(normalise_2d_(frequencies_2d_(x, 0, 1))))
  expect_equal(mutinf10(x, 1, 2), mutinf10(x, 2, 1))
  
  ###########
})

test_that("entropy", {
  expect_equal(entropyE(mtcars[, 1]), entropyE(mtcars[, 1, drop = FALSE]))
  expect_equal(entropyE(mtcars[, 4]), entropyE(mtcars[, 4, drop = FALSE]))
  
  expect_equal(entropy2(mtcars[, 1]), entropy2(mtcars[, 1, drop = FALSE]))
  expect_equal(entropy2(mtcars[, 4]), entropy2(mtcars[, 4, drop = FALSE]))
  
  expect_equal(entropy10(mtcars[, 1]), entropy10(mtcars[, 1, drop = FALSE]))
  expect_equal(entropy10(mtcars[, 4]), entropy10(mtcars[, 4, drop = FALSE]))
})

test_that("Conditional entropy", {
  expect_equal(entropy_condE(mtcars, 1, 4), 
               entropyE(mtcars[, c(1, 4)]) - entropyE(mtcars[, 4, FALSE]))
  expect_equal(entropy_condE(mtcars, 4, 1), 
               entropyE(mtcars[, c(1, 4)]) - entropyE(mtcars[, 1, FALSE]))
  
  
  expect_equal(entropy_cond2(mtcars, 1, 4), 
               entropy2(mtcars[, c(1, 4)]) - entropy2(mtcars[, 4, FALSE]))
  expect_equal(entropy_cond2(mtcars, 4, 1), 
               entropy2(mtcars[, c(1, 4)]) - entropy2(mtcars[, 1, FALSE]))
  
  
  expect_equal(entropy_cond10(mtcars, 1, 4), 
               entropy10(mtcars[, c(1, 4)]) - entropy10(mtcars[, 4, FALSE]))
  expect_equal(entropy_cond10(mtcars, 4, 1), 
               entropy10(mtcars[, c(1, 4)]) - entropy10(mtcars[, 1, FALSE]))
})



test_that("infdist", {
  idx_x <- match("mpg", colnames(mtcars))
  idx_y <- match(c("hp", "wt"), colnames(mtcars))
  
  
  # E
  H_joint <- entropyE(mtcars[, c("mpg", "hp", "wt")])
  H_x <- entropyE(mtcars[, "mpg"])
  H_y <- entropyE(mtcars[, c("hp", "wt")])
  I_xy <- mutinfE(mtcars, idx_x, idx_y)
  H_x_y <- H_x - I_xy
  expect_equal(H_x_y, entropy_condE(mtcars, idx_x, idx_y))
  H_y_x <- H_y - I_xy
  expect_equal(H_y_x, entropy_condE(mtcars, idx_y, idx_x))
  id <- (H_x_y + H_y_x) / H_joint
  expect_equal(id, infdistE(mtcars, idx_x, idx_y))
  expect_equal(id, infdistE(mtcars, idx_y, idx_x))
  expect_equal(id, 0.1)
  
  
  # 2
  H_joint <- entropy2(mtcars[, c("mpg", "hp", "wt")])
  H_x <- entropy2(mtcars[, "mpg"])
  H_y <- entropy2(mtcars[, c("hp", "wt")])
  I_xy <- mutinf2(mtcars, idx_x, idx_y)
  H_x_y <- H_x - I_xy
  expect_equal(H_x_y, entropy_cond2(mtcars, idx_x, idx_y))
  H_y_x <- H_y - I_xy
  expect_equal(H_y_x, entropy_cond2(mtcars, idx_y, idx_x))
  id <- (H_x_y + H_y_x) / H_joint
  expect_equal(id, infdist2(mtcars, idx_x, idx_y))
  expect_equal(id, infdist2(mtcars, idx_y, idx_x))
  expect_equal(id, 0.1)
  
  # 10
  H_joint <- entropy10(mtcars[, c("mpg", "hp", "wt")])
  H_x <- entropy10(mtcars[, "mpg"])
  H_y <- entropy10(mtcars[, c("hp", "wt")])
  I_xy <- mutinf10(mtcars, idx_x, idx_y)
  H_x_y <- H_x - I_xy
  expect_equal(H_x_y, entropy_cond10(mtcars, idx_x, idx_y))
  H_y_x <- H_y - I_xy
  expect_equal(H_y_x, entropy_cond10(mtcars, idx_y, idx_x))
  id <- (H_x_y + H_y_x) / H_joint
  expect_equal(id, infdist10(mtcars, idx_x, idx_y))
  expect_equal(id, infdist10(mtcars, idx_y, idx_x))
  expect_equal(id, 0.1)
})
