test_that("vcf_to_gds works", {
  expect_error(vcf_to_gds(3))
})

test_that("gds_to_pedmap works", {
  expect_error(gds_to_pedmap(3))
})
