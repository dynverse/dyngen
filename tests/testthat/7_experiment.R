context("7_experiment")

network <- tribble(
  ~from, ~to, ~length, ~directed,
  "a", "b", 1, TRUE,
  "b", "c", 2, TRUE,
  "c", "d", 1, TRUE
)
sim_meta <- tribble(
  ~simulation_i, ~sim_time, ~from, ~to, ~time, ~orig_ix,
  1, 0, "a", "b", 0, 1,
  1, .1, "a", "b", .25, 2,
  1, .2, "a", "b", .5, 3,
  1, .3, "a", "b", .75, 4,
  1, .4, "a", "b", 1, 5,
  1, .5, "b", "c", 0, 6,
  1, .6, "b", "c", 0.25, 7,
  1, .7, "b", "c", 0.5, 8,
  1, .8, "b", "c", 0.75, 9,
  1, .9, "b", "c", 1, 10
)
num_cells <- 5

test_that("experiment_synchronised divides cells up correctly", {
  params <- experiment_synchronised()
  
  ix <- .generate_experiment_synchronised(
    network = network,
    sim_meta = sim_meta,
    params = params,
    num_cells = num_cells
  )
  
  expect_equal(length(ix), num_cells)
})


test_that("experiment_snapshot divides cells up correctly", {
  params <- experiment_snapshot()
  
  expect_warning({
    ix <- .generate_experiment_snapshot(
      network = network,
      sim_meta = sim_meta,
      params = params,
      num_cells = num_cells
    )
  })
  
  expect_equal(length(ix), num_cells)
})