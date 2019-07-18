context(".generate_partitions")

test_that("num_tfs_sampler works as expected", {
  num_tfs <- 90
  num_modules <- 8
  
  samples <- map(
    seq_len(10000),
    ~ .generate_partitions(num_tfs, num_modules, min_elements_per_group = 1)
  )
  
  expect_true(all(map_int(samples, length) == num_modules))
  expect_true(all(map_int(samples, sum) == num_tfs))
  expect_true(all(map_int(samples, min) >= 1))
  
  avg <- Reduce("+", samples) / length(samples)
  
  exp <- num_tfs / num_modules
  expect_true(all(abs(avg - exp) / exp < 0.1))
})

test_that("num_tfs_sampler works as expected", {
  num_tfs <- 200
  num_modules <- 13
  
  samples <- map(
    seq_len(10000),
    ~ .generate_partitions(num_tfs, num_modules, min_elements_per_group = 1)
  )
  
  expect_true(all(map_int(samples, length) == num_modules))
  expect_true(all(map_int(samples, sum) == num_tfs))
  expect_true(all(map_int(samples, min) >= 1))
  
  avg <- Reduce("+", samples) / length(samples)
  
  exp <- num_tfs / num_modules
  expect_true(all(abs(avg - exp) / exp < 0.1))
})



test_that("num_tfs_sampler works as expected", {
  num_tfs <- 10000
  num_modules <- 8
  
  samples <- map(
    seq_len(10000),
    ~ .generate_partitions(num_tfs, num_modules, min_elements_per_group = 1)
  )
  
  expect_true(all(map_int(samples, length) == num_modules))
  expect_true(all(map_int(samples, sum) == num_tfs))
  expect_true(all(map_int(samples, min) >= 1))
  
  avg <- Reduce("+", samples) / length(samples)
  
  exp <- num_tfs / num_modules
  expect_true(all(abs(avg - exp) / exp < 0.1))
})


test_that("num_tfs_sampler works as expected", {
  expect_equal(.generate_partitions(50, 10, 5), rep(5, 10))
})


