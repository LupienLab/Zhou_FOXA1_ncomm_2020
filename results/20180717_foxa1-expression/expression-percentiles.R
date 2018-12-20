# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read TCGA data from the Xena RDS file
tcga_xena = readRDS("../../data/external/TCGA-PRAD/TCGA_Xena_PRAD.rds")

# read DEPMAP expression data
depmap = fread(
    "../../data/external/DEPMAP/CCLE_depMap_18Q4_TPM_v2.csv",
    header = TRUE
)
colnames(depmap)[1] = "Line"
# read DEPMAP metadata
metadata = fread(
    "../../data/external/DEPMAP/DepMap-2018q4-celllines.csv",
    header = TRUE
)

# select prostate lines
prostate_lines = metadata[
    get("Primary Disease") == "Prostate Cancer",
    .(DepMap_ID, Aliases)
]
colnames(prostate_lines) = c("ID", "Name")
# map cell names to more visually-friendly names, consistent with other plots
prostate_lines[1, Name := "PC3"]
prostate_lines[5, Name := "MDA PCa 2B"]
prostate_lines[6, Name := "22Rv1"]
prostate_lines[7, Name := "LNCaP"]
prostate_lines[8, Name := "DU145"]

pca_exprs = melt(
    depmap[Line %in% prostate_lines[, ID], ],
    id.vars = 1,
    variable.name = "Gene",
    value.name = "log2TPM"
)
# add cell line names
pca_exprs = merge(
    pca_exprs,
    prostate_lines,
    by.x = "Line",
    by.y = "ID"
)

# get FOXA1 expression in prostate lines
foxa1_lines = pca_exprs[grep("FOXA1", Gene), .(Name, log2TPM)]

# ==============================================================================
# Analysis
# ==============================================================================
# find FOXA1 expression percentile wrt all other mapped reads for that patient
percentiles_tcga = apply(tcga_xena, 2, function(samp) {
    trunc(rank(samp, na.last = NA)) / sum(!is.na(samp))
})

percentiles_df = as.data.frame(percentiles_tcga)

# save data
fwrite(
    percentiles_df,
    "percentiles-tcga.tsv",
    col.names = TRUE,
    sep = "\t",
    row.names = TRUE
)

# find FOXA1 expression percentiles in DEPMAP cell line data
percentiles = sapply(
    pca_exprs[, unique(Name)],
    function(n) {
        percentile = ecdf(pca_exprs[Name == n, log2TPM])
        percentile(foxa1_exprs[Name == n, log2TPM])
    }
)

# save percentiles as table
dt = data.table(
    Line = names(percentiles),
    Percentile = percentiles
)
fwrite(
    dt,
    "percentiles-lines.tsv",
    col.names = TRUE,
    sep = "\t"
)
