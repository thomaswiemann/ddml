test_that("parse_parallel handles NULL inputs safely", {
  res <- parse_parallel(NULL)
  expect_equal(res$num_cores, 1)
  expect_null(res$export)
  expect_null(res$packages)
})

test_that("parse_parallel extracts list structures correctly", {
  pl <- list(cores = 4, export = c("var1"), packages = c("utils"))
  res <- parse_parallel(pl)
  
  expect_equal(res$num_cores, 4)
  expect_equal(res$export, "var1")
  expect_equal(res$packages, "utils")
})

test_that("parse_parallel catches malformed arguments", {
  expect_error(parse_parallel(4), "'parallel' must be a list or NULL")
  expect_error(parse_parallel(TRUE), "'parallel' must be a list or NULL")
})
