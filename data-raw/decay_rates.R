library(tidyverse)

# download data from Schwannhäusser et al., 2011, doi.org/10.1038/nature10098
file <- "data-raw/input_rates_schwannhausser2011.tsv"
if (!file.exists(file)) {
  tmpfile <- tempfile()
  download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fnature10098/MediaObjects/41586_2011_BFnature10098_MOESM304_ESM.xls", tmpfile)
  tab <- readxl::read_excel(tmpfile)
  
  tab2 <- tab %>% 
    select(
      protein_ids = `Protein IDs`,
      protein_names = `Protein Names`,
      gene_names = `Gene Names`,
      uniprot_ids = `Uniprot IDs`,
      refset_protein_ids = `RefSeq protein IDs`,
      refseq_mrna_id = `Refseq mRNA ID`,
      ensembl_id = `ENSEMBL ID`,
      mgi_id = `MGI ID`,
      protein_length = `Protein length [amino acids]`,
      protein_copynumber = `Protein copy number average [molecules/cell]`,
      protein_halflife = `Protein half-life average [h]`,
      mrna_copynumber = `mRNA copy number average [molecules/cell]`,
      mrna_halflife = `mRNA half-life average [h]`,
      transcription_rate = `transcription rate (vsr) average [molecules/(cell*h)]`,
      translation_rate = `translation rate constant (ksp) average [molecules/(mRNA*h)]`
    )
  write_tsv(tab2, file)
}

tab <- read_tsv(file)

# perform data imputation
out <- mice::mice(tab %>% select_if(is.numeric), m = 1, maxit = 10, method = "pmm", seed = 1)

mmu_imputed_rates <- mice::complete(out)

usethis::use_data(mmu_imputed_rates)


# mRNA half lives ---------------------------------------------------------
mrna_halflife_hours <- tab$mrna_halflife %>% na.omit() %>% as.numeric()

dis <- fitdistrplus::fitdist(mrna_halflife_hours, "lnorm")

sam <- rlnorm(10000, meanlog = dis$estimate[["meanlog"]], sdlog = dis$estimate[["sdlog"]])
his <- hist(mrna_halflife_hours, breaks = 100)

den <- density(sam)
den$y <- den$y / max(den$y) * max(his$counts)
lines(den, col = "blue")

dis$estimate


# protein half lives ------------------------------------------------------
protein_halflife_hours <- tab$protein_halflife %>% na.omit() %>% as.numeric()
dis <- fitdistrplus::fitdist(protein_halflife_hours, "lnorm")

sam <- rlnorm(10000, meanlog = dis$estimate[["meanlog"]], sdlog = dis$estimate[["sdlog"]])
his <- hist(protein_halflife_hours, breaks = 100)

den <- density(sam)
den$y <- den$y / max(den$y) * max(his$counts)
lines(den, col = "blue")

dis$estimate


# transcription rate ------------------------------------------------------
transcription_rate <- tab$transcription_rate %>% na.omit() %>% as.numeric()
dis <- fitdistrplus::fitdist(transcription_rate, "lnorm")

sam <- rlnorm(10000, meanlog = dis$estimate[["meanlog"]], sdlog = dis$estimate[["sdlog"]])
his <- hist(transcription_rate, breaks = 100)

den <- density(sam)
den$y <- den$y / max(den$y) * max(his$counts)
lines(den, col = "blue")

dis$estimate

# transcription rate ------------------------------------------------------
translation_rate <- tab$translation_rate %>% na.omit() %>% as.numeric()
dis <- fitdistrplus::fitdist(translation_rate, "lnorm")

sam <- rlnorm(10000, meanlog = dis$estimate[["meanlog"]], sdlog = dis$estimate[["sdlog"]])
his <- hist(translation_rate, breaks = 100)

den <- density(sam)
den$y <- den$y / max(den$y) * max(his$counts)
lines(den, col = "blue")

dis$estimate






# Ratio between mRNA and proteins -----------------------------------------

# The ratio between mRNA and proteins should be around 0.1-0.5
# https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-4-30/tables/1
