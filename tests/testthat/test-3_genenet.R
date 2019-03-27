context("3_genenet")

test_that("num_tfs_sampler works as expected", {
  num_tfs <- 90
  num_modules <- 8
  
  samples <- map(
    seq_len(10000),
    ~ num_tfs_sampler(num_tfs = num_tfs, num_modules = num_modules)
  )
  
  expect_true(all(map_int(samples, length) == num_modules))
  expect_true(all(map_int(samples, sum) == num_tfs))
  expect_true(all(map_int(samples, min) >= 1))
  
  avg <- Reduce("+", samples) / length(samples)
  
  expect_true(all(abs(avg - num_tfs / num_modules) < 1))
})

test_that("num_tfs_sampler works as expected", {
  num_tfs <- 200
  num_modules <- 13
  
  samples <- map(
    seq_len(10000),
    ~ num_tfs_sampler(num_tfs = num_tfs, num_modules = num_modules)
  )
  
  expect_true(all(map_int(samples, length) == num_modules))
  expect_true(all(map_int(samples, sum) == num_tfs))
  expect_true(all(map_int(samples, min) >= 1))
  
  avg <- Reduce("+", samples) / length(samples)
  
  expect_true(all(abs(avg - num_tfs / num_modules) < 1))
})
