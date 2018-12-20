# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
rnai = fread(
    "essentiality-scores.tsv",
    header = TRUE,
    sep = "\t"
)

prostate = rnai[Tissue == "Prostate"]
foxa1_prostate = rnai[Tissue == "Prostate" & Gene == "FOXA1"]

# ==============================================================================
# Analysis
# ==============================================================================
percentiles = sapply(
    prostate[, unique(Line)],
    function(l) {
        percentile = ecdf(prostate[Line == l, Score])
        percentile(foxa1_prostate[Line == l, Score])
    }
)

dt = data.table(
    Line = names(percentiles),
    Percentile = percentiles
)
fwrite(
    dt,
    "percentiles-foxa1.tsv",
    col.names = TRUE,
    sep = "\t"
)
