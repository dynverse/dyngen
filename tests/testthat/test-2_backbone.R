context("2_backbone")

test_that("Testing normal use case of backbone creation", {
  backbone <- backbone(
    module_info = tibble(
      module_id = c("A", "B"),
      basal = c(0.2, 0.8),
      burn = c(TRUE, FALSE),
      independence = c(0.4, 0.5)
    ),
    module_network = tibble(
      from = c("A", "A"),
      to = c("A", "B"),
      effect = c(-1L, 1L),
      strength = c(10, 1),
      cooperativity = c(2, 3)
    ),
    expression_patterns = tibble(
      from = c("burn", "start"),
      to = c("start", "end"),
      module_progression = c("+A", "+B"),
      start = c(TRUE, FALSE),
      burn = c(TRUE, FALSE),
      time = c(1, 1)
    )
  )
  
  expect_is(backbone, "list")
  expect_named(backbone, c("module_info", "module_network", "expression_patterns"))
  
  expect_is(backbone$module_info, "data.frame")
  expect_named(backbone$module_info, c("module_id", "basal", "burn", "independence", "color"))
  expect_equal(backbone$module_info$module_id, c("A", "B"))
  expect_equal(backbone$module_info$basal, c(.2, .8))
  expect_equal(backbone$module_info$burn, c(TRUE, FALSE))
  expect_equal(backbone$module_info$independence, c(.4, .5))
  expect_is(backbone$module_info$color, "character")
  
  expect_is(backbone$module_network, "data.frame")
  expect_equal(backbone$module_network$from, c("A", "A"))
  expect_equal(backbone$module_network$to, c("A", "B"))
  expect_equal(backbone$module_network$effect, c(-1L, 1L))
  expect_equal(backbone$module_network$strength, c(10, 1))
  expect_equal(backbone$module_network$cooperativity, c(2, 3))
  
  expect_is(backbone$expression_patterns, "data.frame")
  expect_equal(backbone$expression_patterns$from, c("burn", "start"))
  expect_equal(backbone$expression_patterns$to, c("start", "end"))
  expect_equal(backbone$expression_patterns$module_progression, c("+A,+B", "+B"))
  expect_equal(backbone$expression_patterns$start, c(TRUE, FALSE))
  expect_equal(backbone$expression_patterns$burn, c(TRUE, FALSE))
  expect_equal(backbone$expression_patterns$time, c(1, 1))
})

# TODO: add test for when things go wrong