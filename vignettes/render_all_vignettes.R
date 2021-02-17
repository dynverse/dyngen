library(tidyverse)

rmds <- list.files("vignettes", pattern = "\\.Rmd", full.names = TRUE)

for (rmd in rmds) {
  cat("Rendering ", rmd, "\n", sep = "")
  start <- Sys.time()
  rmarkdown::render(rmd, output_format = rmarkdown::github_document(), quiet = TRUE)
  end <- Sys.time()
  cat("  Time elapsed: ", as.numeric(difftime(end, start, units = "secs")), "s\n", sep = "")
}
