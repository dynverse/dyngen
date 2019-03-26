generator_config <- function(
  modulenet = modulenet_linear(),
  platform = platform_simple(),
  targetadder = targetadder_realnet()
) {
  lst(
    modulenet,
    platform,
    targetadder
  )
}