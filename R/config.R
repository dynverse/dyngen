generator_config <- function(
  modulenet = modulenet_linear(),
  platform = platform_simple(),
  tf_generator = tf_random(),
  target_generator = target_realnet(),
  kinetics_generator = kinetics_samplers()
) {
  lst(
    modulenet,
    platform,
    tf_generator,
    target_generator,
    kinetics_generator
  )
}