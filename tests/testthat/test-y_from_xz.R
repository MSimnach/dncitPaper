test_that("basic_equalities", {
  s <- c(1,2,3,4,5)

  g <- y_from_xz(post_non_lin = 1)
  expect_equal(g(s), s)

  g <- y_from_xz(post_non_lin = 2)
  expect_equal(g(s), s^2)

  g <- y_from_xz(post_non_lin = 3)
  expect_equal(g(s), s^3)

  t <- scale(s)
  g <- y_from_xz(post_non_lin = 4)
  expect_equal(g(s), exp(-t^2/2) * sin(3 * t))

  g <- y_from_xz(post_non_lin = 5)
  expect_equal(g(s), exp(-t^2/2) * sin(24 * t))
})

