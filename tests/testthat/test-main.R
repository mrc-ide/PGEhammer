test_that("square() handles default/null values", {
  expect_equal(square(), square(1:5))
  expect_error(square(NULL))
})
