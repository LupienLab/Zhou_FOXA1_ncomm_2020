# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read TCGA manifest data
tcga <- fread(
    "tcga-aggregated.tsv",
    header = TRUE,
    sep = "\t"
)

# load Ensembl gene IDs
ensembl <- fread(
    "ensembl-gene-IDs.tsv",
    header = TRUE,
    sep = "\t",
    select = 1:2
)

# map Ensembl ID to a gene name
tcga <- merge(
    x = tcga,
    y = ensembl,
    by.x = "EnsemblID",
    by.y = "ensembl_gene_id",
    all.x = TRUE
)

# ==============================================================================
# Analysis
# ==============================================================================
# find FOXA1 expression percentile wrt all other mapped reads for that patient
tcga_percentile <- as.data.table(apply(
    tcga[, c(-1, -552)],  # only operate on patient data
    2,
    function(x) {
        ifelse(
            is.na(x),
            NA,
            trunc(rank(x, na.last = NA))/sum(!is.na(x))
        )
    }
))
tcga_percentile[, EnsemblID := tcga$EnsemblID]
tcga_percentile[, Description := tcga$hgnc_symbol]
# save data
fwrite(
    tcga_percentile,
    "tcga-expression-percentiles.tsv",
    col.names = TRUE,
    sep = "\t"
)