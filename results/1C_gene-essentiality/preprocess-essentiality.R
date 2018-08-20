# ==============================================================================
# Environment
# ==============================================================================
suppressWarnings(library("data.table"))
suppressWarnings(library("argparse"))
suppressWarnings(library("stringr"))

if (!interactive()) {
    PARSER <- argparse::ArgumentParser(
        description = "Preprocess Achilles data for simpler processing"
    )
    PARSER$add_argument(
        "rnai",
        type = "character",
        help = "GCT file containing Achilles RNAi-based essentiality scores"
    )
    PARSER$add_argument(
        "crispr",
        type = "character",
        help = "GCT file containing Achilles CRISPR-based essentiality scores"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        rnai = "../../data/external/Achilles_QC_v2.4.3.rnai.Gs.gct",
        crispr = "../../data/external/Achilles_v3.3.8.Gs.gct"
    )
}

# ==============================================================================
# Data
# ==============================================================================
achilles_crispr <- fread(
    ARGS$crispr,
    header = TRUE,
    skip = 2
)
achilles_rnai <- fread(
    ARGS$rnai,
    header = TRUE,
    skip = 2
)

# Cell lines of interest
crispr_lines <- grep("PROSTATE", colnames(achilles_crispr))
rnai_lines <- grep("PROSTATE", colnames(achilles_rnai))

prostate_crispr <- c(
    "LNCAPCLONEFGC_PROSTATE",
    "PC3_PROSTATE"
)
prostate_rnai <- c(
    "22RV1_PROSTATE",
    "NCHIH660_PROSTATE",
    "VCAP_PROSTATE"
)
gene <- "FOXA1"

# melted data.table for essentiality across all cell lines
crispr_table <- melt(
    achilles_crispr,
    id.vars = "Description",
    measure.vars = 3:35,
    variable.name = "Cell",
    value.name = "Score"
)
# find index of underscore in string
crispr_table[, Underscore := as.vector(regexpr("_", Cell))]
# classify cell line based on tissue type
crispr_table[, Tissue := substring(Cell, Underscore + 1, 100)]
crispr_table[, Tissue := factor(str_to_title(Tissue))]
# extract readable name of cell line
crispr_table[, Line := substring(Cell, 1, Underscore - 1)]
crispr_table[, Line := factor(Line)]
# remove unused column
crispr_table[, Underscore := NULL]

rnai_table <- melt(
    achilles_rnai,
    id.vars = "Description",
    measure.vars = 3:218,
    variable.name = "Cell",
    value.name = "Score"
)
# find index of underscore in string
rnai_table[, Underscore := as.vector(regexpr("_", Cell))]
# classify cell line based on tissue type
rnai_table[, Tissue := substring(Cell, Underscore + 1, 100)]
rnai_table[, Tissue := factor(str_to_title(Tissue))]
# extract readable name of cell line
rnai_table[, Line := substring(Cell, 1, Underscore - 1)]
rnai_table[, Line := factor(Line)]
# remove unused column
rnai_table[, Underscore := NULL]

# ==============================================================================
# Save
# ==============================================================================
fwrite(
    rnai_table,
    "essentiality-rnai.tsv",
    col.names = TRUE,
    sep = "\t"
)
fwrite(
    crispr_table,
    "essentiality-crispr.tsv",
    col.names = TRUE,
    sep = "\t"
)
