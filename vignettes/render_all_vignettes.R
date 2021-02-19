library(tidyverse)

articles <- tibble(
  rmd = list.files("vignettes", pattern = "\\.Rmd", full.names = TRUE, recursive = TRUE),
  md = gsub("Rmd$", "md", rmd),
  rmd_mtime = file.info(rmd)$mtime,
  md_mtime = file.info(md)$mtime,
  render = is.na(md_mtime) | rmd_mtime > md_mtime,
  html_preview = !grepl("advanced_", rmd) && !grepl("showcase_backbones", rmd)
)

pwalk(articles %>% filter(render), function(rmd, md, html_preview, ...) {
  cat("Rendering ", rmd, "\n", sep = "")
  start <- Sys.time()
  rmarkdown::render(rmd, output_format = rmarkdown::github_document(html_preview = html_preview), quiet = TRUE)
  end <- Sys.time()
  cat("  Time elapsed: ", as.numeric(difftime(end, start, units = "secs")), "s\n", sep = "")
})
