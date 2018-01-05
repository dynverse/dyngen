context("Test random tree generation")

test_that("A random tree can be generated and converted to a module network", {
  stagenet <- generate_random_tree()
  modulenet <- from_stages_to_modulenet(stagenet)
})