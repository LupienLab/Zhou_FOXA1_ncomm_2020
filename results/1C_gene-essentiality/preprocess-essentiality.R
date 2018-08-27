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
        rnai = "../../data/external/Achilles/Achilles_QC_v2.4.3.rnai.Gs.gct",
        crispr = "../../data/external/Achilles/Achilles_v3.3.8.Gs.gct"
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
# further underscores in Tissue replaced with spaces
crispr_table[, Tissue := gsub("_", " ", Tissue)]
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
# further underscores in Tissue replaced with spaces
rnai_table[, Tissue := gsub("_", " ", Tissue)]
# remove unused column
rnai_table[, Underscore := NULL]

# combine tables for single data file
crispr_table[, Method := "CRISPR"]
rnai_table[, Method := "RNAi"]
all_table <- rbindlist(
    list(crispr_table, rnai_table)
)

# fix the naming on a few selected cell lines
all_table[Line == "22RV1", Line := "22Rv1"]
all_table[Line == "LNCAPCLONEFGC", Line := "LNCaP"]
all_table[Line == "NCIH660", Line := "NCI-H660"]
all_table[Line == "VCAP", Line := "VCaP"]

# ==============================================================================
# Save
# ==============================================================================
fwrite(
    all_table,
    "essentiality-scores.tsv",
    col.names = TRUE,
    sep = "\t"
)
