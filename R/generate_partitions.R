.generate_partitions <- function(num_elements, num_groups, min_elements_per_group) {
  assert_that(
    min_elements_per_group >= 0,
    num_groups * min_elements_per_group <= num_elements
  )

  free <- num_elements - num_groups * min_elements_per_group
  sample(seq(0, free), num_groups - 1, replace = TRUE) |>
    sort() |>
    (\(x) c(0, x, free))() |>
    diff() |>
    (\(x) x + min_elements_per_group)() |>
    as.integer()
}
