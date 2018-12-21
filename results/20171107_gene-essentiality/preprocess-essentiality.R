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
        "essentiality",
        type = "character",
        help = "CSV file containing DEPMAP gene essentiality scores"
    )
    PARSER$add_argument(
        "metadata",
        type = "character",
        help = "Metadata for cell lines used"
    )
    ARGS <- PARSER$parse_args()
} else {
    ARGS <- list(
        essentiality = "../../data/external/DEPMAP/D2_combined_gene_dep_scores.csv",
        metadata = "../../data/external/DEPMAP/DepMap-2018q4-celllines.csv"
    )
}

# ==============================================================================
# Data
# ==============================================================================
# read in RNAi data
cat("Reading data\n")
rnai = fread(ARGS$essentiality, header = TRUE)
# change first column name
colnames(rnai)[1] = "Gene"

# read in cell line metadata
metadata = fread(ARGS$metadata, header = TRUE)

# ==============================================================================
# Preprocessing
# ==============================================================================
cat("Removing non-cancer cell lines\n")
# filter out non-cancer cell lines
#   list of all non-cancer values in "Primary Disease" column in metadata
noncancer_disease_classes = c(
    "Fibroblast",
    "Immortalized",
    "immortalized_epithelial",
    "Non-Cancerous",
    "Primary Cells",
    "unknown"
)
#   find all non-cancerous lines
noncancer_lines = metadata[
    get("Primary Disease") %in% noncancer_disease_classes,
    CCLE_Name
]
#   remove columns corresponding to all non-cancerous lines
noncancer_col_idx = which(colnames(rnai) %in% noncancer_lines)
rnai = rnai[, .SD, .SDcols = -noncancer_col_idx]

cat("Aggregating and parsing data\n")
rnai_table = melt(
    rnai,
    id.vars = "Gene",
    variable.name = "Cell",
    value.name = "Score"
)
# clean up text for simpler parsing
#   find index of underscore in string
rnai_table[, Underscore := as.vector(regexpr("_", Cell))]
#   classify cell line based on tissue type
rnai_table[, Tissue := substring(Cell, Underscore + 1, 100)]
rnai_table[, Tissue := factor(str_to_title(Tissue))]
#   extract readable name of cell line
rnai_table[, Line := substring(Cell, 1, Underscore - 1)]
rnai_table[, Line := factor(Line)]
#   further underscores in Tissue replaced with spaces
rnai_table[, Tissue := gsub("_", " ", Tissue)]
#   remove unused column
rnai_table[, Underscore := NULL]
#   find location of brackets in Gene column
rnai_table[, Gene := gsub(" \\(\\d+\\)", "", Gene)]

# map cell names to more visually-friendly names, consistent with other plots
rnai_table[Line == "22RV1", Line := "22Rv1"]
rnai_table[Line == "LNCAPCLONEFGC", Line := "LNCaP"]
rnai_table[Line == "MDAPCA2B", Line := "MDA PCa 2B"]
rnai_table[Line == "NCIH660", Line := "NCI-H660"]
rnai_table[Line == "PRECLH", Line := "PrEC LH"]
rnai_table[Line == "VCAP", Line := "VCaP"]

# ==============================================================================
# Save
# ==============================================================================
cat("Saving combined data\n")
fwrite(
    rnai_table,
    "essentiality-scores.tsv",
    col.names = TRUE,
    sep = "\t"
)
