library(tidyverse)

articles <- tibble(
  rmd = list.files("vignettes", pattern = "\\.Rmd", full.names = TRUE),
  md = gsub("Rmd$", "md", rmd),
  rmd_mtime = file.info(rmd)$mtime,
  md_mtime = file.info(md)$mtime,
  render = is.na(md_mtime) | rmd_mtime > md_mtime
)

pwalk(articles %>% filter(render), function(rmd, md, ...) {
  cat("Rendering ", rmd, "\n", sep = "")
  start <- Sys.time()
  rmarkdown::render(rmd, output_format = rmarkdown::github_document(html_preview = FALSE), quiet = TRUE)
  end <- Sys.time()
  cat("  Time elapsed: ", as.numeric(difftime(end, start, units = "secs")), "s\n", sep = "")
})
