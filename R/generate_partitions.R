.generate_partitions <- function(num_elements, num_groups, min_elements_per_group) {
  # satisfy r cmd check
  `.` <- NULL
  
  assert_that(
    min_elements_per_group >= 0, 
    num_groups * min_elements_per_group <= num_elements
  )
  
  sample(
    seq(0, num_elements - num_groups * min_elements_per_group),
    num_groups - 1,
    replace = TRUE
  ) %>% 
    sort() %>% 
    c(0, ., num_elements - num_groups * min_elements_per_group) %>% 
    diff() %>% 
    { . + min_elements_per_group} %>% 
    as.integer()
}