
modnet_edge <- function(from, to, effect = 1L, strength = 1, hill = 2) {
  max_len <- max(length(from), length(to))
  assert_that(
    length(from) == 1 || length(to) == 1 || length(from) == length(to),
    length(effect) == 1 || length(effect) == max_len,
    length(strength) == 1 || length(strength) == max_len,
    length(hill) == 1 || length(hill) == max_len
  )
  tibble(
    from = from,
    to = to,
    effect = effect,
    strength = strength,
    hill = hill
  )
}
modnet_chain <- function(mids, effect = 1L, strength = 1, hill = 2) {
  if (length(mids) == 1) return(NULL)
  modnet_edge(
    from = mids %>% head(-1),
    to = mids %>% tail(-1),
    effect = effect,
    strength = strength,
    hill = hill
  )
}
modnet_self <- function(mids, effect = 1L, strength = 1, hill = 2) {
  modnet_edge(
    from = mids,
    to = mids,
    effect = effect,
    strength = strength,
    hill = hill
  )
}
modnet_pairwise <- function(from, to = from, self = FALSE, effect = 1L, strength = 1, hill = 2) {
  # satisfy r cmd check
  fromi <- toi <- NULL

  max_len <- length(from) * ifelse(self, length(to), length(to) - 1)
  assert_that(
    self || length(from) == length(to),
    length(effect) == 1 || length(effect) == max_len,
    length(strength) == 1 || length(strength) == max_len,
    length(hill) == 1 || length(hill) == max_len
  )
  df <- crossing(
    fromi = seq_along(from),
    toi = seq_along(to)
  ) 
  if (!self) {
    df <- df %>% filter(fromi != toi)
  }
  df %>% 
    transmute(
      from = from[fromi],
      to = to[toi],
      effect,
      strength,
      hill
    )
}