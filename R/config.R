#' @export
generator_config <- function(
  modulenet = modulenet_linear(),
  platform = platform_simple(),
  tfgen_params = tf_random(),
  targetgen_params = target_realnet(),
  kinetics_params = kinetics_default()
) {
  lst(
    modulenet,
    platform,
    tfgen_params,
    targetgen_params,
    kinetics_params
  )
}