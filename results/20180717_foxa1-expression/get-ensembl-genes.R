suppressMessages(library("biomaRt"))
suppressMessages(library("data.table"))

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
ensembl_ids <- getBM(
    attributes = c(
        'ensembl_gene_id',
        'hgnc_symbol',
        'chromosome_name',
        'start_position',
        'end_position',
        'band'
    ),
    mart = ensembl
)
ensembl_ids <- as.data.table(ensembl_ids)
fwrite(
    x = ensembl_ids,
    file = "ensembl-gene-IDs.tsv",
    sep = "\t",
    col.names = TRUE
)
