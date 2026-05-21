# =========================================================
# Script to calculate methylation accuracy metrics
# using ONT methylation calls and bisulfite sequencing data
#
# Usage:
# Rscript methylation_metrics.R <ont_file.tsv> <bis_file.tsv> <motif>
#
# Example:
# Rscript methylation_metrics.R \
#   ont_calls.tsv \
#   bisulfite.tsv \
#   CG
# =========================================================

# -------------------------
# Load required libraries
# -------------------------
library(data.table,warn.conflicts = F)
library(dplyr,warn.conflicts = F)

# -------------------------
# Read command line arguments
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments are provided
if(length(args) != 3){
  stop(
    paste0(
      "Usage: Rscript methylation_metrics.R ",
      "<ont_file.tsv> <bis_file.tsv> <motif>"
    )
  )
}

# Assign command line arguments
ont_file <- args[1]
bis_file <- args[2]
motif <- args[3]

# -------------------------
# Read ONT methylation calls
# -------------------------
ont_df <- fread(
  ont_file,
  sep = '\t',
  header = TRUE,
  stringsAsFactors = TRUE
)

# -------------------------
# Read bisulfite sequencing file
# -------------------------
bis_df <- fread(
  bis_file,
  sep = '\t',
  header = FALSE,
  stringsAsFactors = TRUE,
  col.names = c(
    "chrom",
    "p1",
    "p2",
    "M",
    "UM",
    "strand",
    "coverage",
    "per",
    "context",
    "fullcontext"
  )
)

# -------------------------
# Filter ONT calls by coverage
# -------------------------
ont_df <- ont_df %>%
  filter(coverage >= 20)

# -------------------------
# Create fully methylated sites
# from bisulfite data
# -------------------------
bis_df_per100 <- bis_df %>%
  filter(
    coverage >= 20,
    per == 100,
    context == motif
  )

# -------------------------
# Create fully unmethylated sites
# from bisulfite data
# -------------------------
bis_df_per0 <- bis_df %>%
  filter(
    coverage >= 20,
    per == 0,
    context == motif
  )

# =========================================================
# TRUE POSITIVE / FALSE NEGATIVE calculation
# =========================================================

set.seed(123)

pos_sites <- bis_df_per100 %>%
  inner_join(
    ont_df,
    by = c("chrom", "p1", "p2", "strand")
  ) %>%
  summarise(
    Species = unique(species),
    Tool = unique(tool),
    Context = unique(context),
    Model = unique(model),
    Sample = unique(sample),
    Count = n(),
    TP = sum(M.y),
    FN = sum(UM.y),
    .groups = "drop"
  )

# =========================================================
# FALSE POSITIVE / TRUE NEGATIVE calculation
# =========================================================

set.seed(123)

neg_sites <- bis_df_per0 %>%
  inner_join(
    ont_df,
    by = c("chrom", "p1", "p2", "strand")
  ) %>%
  slice_sample(n = pos_sites$Count) %>%
  summarise(
    Species = unique(species),
    Tool = unique(tool),
    Context = unique(context),
    Model = unique(model),
    Sample = unique(sample),
    Count = n(),
    FP = sum(M.y),
    TN = sum(UM.y),
    .groups = "drop"
  )

# =========================================================
# Calculate performance metrics
# =========================================================

accuracy_metrics <- pos_sites %>%
  inner_join(
    neg_sites,
    by = c(
      "Species",
      "Tool",
      "Count",
      "Context",
      "Model",
      "Sample"
    )
  ) %>%
  mutate(
    Precision = TP / (TP + FP),
    Recall = TP / (TP + FN),
    F1 = (2 * Precision * Recall) / (Precision + Recall),
    Specificity = TN / (TN + FP),
    Accuracy = (TP + TN) / (TP + TN + FP + FN),
    Pos_sites = Count,
    Neg_sites = Count
  ) %>%
  select(
    Species,
    Tool,
    Context,
    Model,
    Sample,
    TP,
    FP,
    TN,
    FN,
    Precision,
    Recall,
    F1,
    Specificity,
    Pos_sites,
    Neg_sites,
    Accuracy
  )

# -------------------------
# Print final metrics
# -------------------------
print(accuracy_metrics,row.names = F)
