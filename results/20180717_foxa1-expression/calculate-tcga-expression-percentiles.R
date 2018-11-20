# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read TCGA data from the Xena RDS file
tcga_xena <- readRDS("../../data/external/TCGA-PRAD/TCGA_Xena_PRAD.rds")

# ==============================================================================
# Analysis
# ==============================================================================
# find FOXA1 expression percentile wrt all other mapped reads for that patient
percentiles = apply(tcga_xena, 2, function(samp) {
    trunc(rank(samp, na.last = NA)) / sum(!is.na(samp))
})

percentiles_df = as.data.frame(percentiles)

# save data
fwrite(
    percentiles_df,
    "tcga-expression-percentiles.tsv",
    col.names = TRUE,
    sep = "\t",
    row.names = TRUE
)
