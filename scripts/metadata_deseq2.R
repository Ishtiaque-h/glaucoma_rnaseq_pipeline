#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

#project_dir <- "glaucoma_rnaseq_pipeline"   # adjust if needed
in_meta     <- file.path("data/metadata", "SRP394552_metadata.csv")
out_meta    <- file.path("data/metadata", "metadata_deseq2.csv")

meta <- read.csv(in_meta, stringsAsFactors = FALSE)

# We keep only the columns we really need for analysis
samples <- meta %>%
  transmute(
    run_id         = run_accession,
    biosample_id   = biosample,
    instrument_model = instrument_model,
    age_raw        = age,
    condition_raw  = condition,  # "IOP" or "control"
    treatment      = treatment
  ) %>%
  # Normalize condition
  mutate(
    condition = tolower(condition_raw)
  )

# Define coarse age group
samples <- samples %>%
  mutate(
    age_group = case_when(
      grepl("^3-month-old \\(at the 1st IOP\\)", age_raw) ~ "3m_rep",  # repeated stress arm (NovaSeq)
      age_raw == "3-month-old"                             ~ "3m",
      age_raw == "18-month-old"                            ~ "18m",
      TRUE ~ NA_character_
    )
  )

# Define stress regime within each age group
samples <- samples %>%
  mutate(
    stress_regime = case_when(
      # Repeated-stress arm (NovaSeq)
      age_group == "3m_rep" & condition == "iop"     ~ "Repeated_IOP",
      age_group == "3m_rep" & condition == "control" ~ "Repeated_Ctrl",

      # Acute vs none (HiSeq arms)
      grepl("unilateral untreated control", treatment)                       ~ "None",
      grepl("constant 30mmHg IOP for 1 hour", treatment)                    ~ "Acute",

      TRUE ~ NA_character_
    )
  )

# Define high-level biological group labels
# We keep repeated controls separate to avoid mixing platforms
samples <- samples %>%
  mutate(
    Group_Label = case_when(
      age_group == "3m"     & stress_regime == "None"         ~ "Young_Ctrl",
      age_group == "3m"     & stress_regime == "Acute"        ~ "Young_Acute",
      age_group == "3m_rep" & stress_regime == "Repeated_IOP" ~ "Young_Rep_IOP",
      age_group == "3m_rep" & stress_regime == "Repeated_Ctrl"~ "Young_Rep_Ctrl",
      age_group == "18m"    & stress_regime == "None"         ~ "Old_Ctrl",
      age_group == "18m"    & stress_regime == "Acute"        ~ "Old_Acute",
      TRUE ~ NA_character_
    )
  )

# Sanity check: every row must have a group
if (any(is.na(samples$Group_Label))) {
  cat("Rows with missing Group_Label:\n")
  print(samples[is.na(samples$Group_Label), ])
  stop("Some samples do not have a Group_Label – fix before proceeding.")
}

# Show counts per group for your inspection
print(table(samples$Group_Label, samples$instrument_model))

# Save simplified metadata
write.csv(samples, out_meta, row.names = FALSE)
cat("Written:", out_meta, "\n")
