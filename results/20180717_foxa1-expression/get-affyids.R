suppressMessages(library("biomaRt"))
suppressMessages(library("data.table"))
# get Affymetrix Probe IDs (RUN WHEN CONNECTED TO INTERNET)

cat("MAKE SURE YOU CAN ACCESS THE INTERNET WHEN RUNNING THIS SCRIPT\n")

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
affyids <- getBM(
    attributes = c(
        'affy_hg_u133_plus_2',
        'hgnc_symbol',
        'chromosome_name',
        'start_position',
        'end_position',
        'band'
    ),
    mart = ensembl
)
affyids <- as.data.table(affyids)
fwrite(
    x = affyids,
    file = "Affymetrix_U133_Plus_2-probe-IDs.tsv",
    sep = "\t",
    col.names = TRUE
)
