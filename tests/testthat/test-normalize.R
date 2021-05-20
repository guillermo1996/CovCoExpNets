test_that("normalize() normalizes a vector", {
  x <- 1:5
  x_norm <- normalize(x)
  y <- (-2:2)/sd(1:5)
  expect_equal(x_norm, c(mean = 3,standard.deviation = sd(1:5),y))
  #expect_equal(x_norm,y)
})
