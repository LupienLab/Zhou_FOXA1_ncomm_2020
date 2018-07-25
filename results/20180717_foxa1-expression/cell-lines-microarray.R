# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("biomaRt"))

# ==============================================================================
# Data
# ==============================================================================
# read CCLE data
ccle <- fread(
    "../../data/external/CCLE/CCLE_Expression_2012-09-29.res",
    header = TRUE,
    sep = "\t",
    blank.lines.skip = TRUE,  # second line of the res file is blank
    fill = TRUE  # third line of res file only has one column
)

# as per the RES format, the 3rd line contains the number of rows in the file
# remove this line from the table
ccle <- ccle[-1]

# as per the RES format, every second column is the Absent/Present/Marginal call
# for the previous column. Ignore these columns and remove them
# (these columns are all preface with 'V' since fread assigns them a dummy name)
cell_line_idx <- grepl("^V\\d+", colnames(ccle))
ccle <- ccle[, which(cell_line_idx) := NULL]

# ==============================================================================
# Analysis
# ==============================================================================
# get Affymetrix Probe IDs (RUN WHEN CONNECTED TO INTERNET)
# (JUST LOAD THE SAVED TABLE WHEN NOT CONNECTED TO INTERNET)
# ensembl <- useMart("ensembl")
# ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
# affyids <- getBM(
#     attributes = c(
#         'affy_hg_u133_plus_2',
#         'hgnc_symbol',
#         'chromosome_name',
#         'start_position',
#         'end_position',
#         'band'
#     ),
#     filters = 'affy_hg_u133_plus_2', 
#     values = ccle$Accession, 
#     mart = ensembl
# )
# affyids <- as.data.table(affyids)
# fwrite(
#     x = affyids,
#     file = "Affymetrix_U133_Plus_2-probe-IDs.tsv",
#     sep = "\t",
#     col.names = TRUE
# )

# load Affymetrix probe IDs and information
affyids <- fread(
    "Affymetrix_U133_Plus_2-probe-IDs.tsv",
    sep = "\t",
    header = TRUE,
    select = 1:2,
    col.names = c("Accession", "Gene")
)

# remove probes that map to multiple annotations as to not skew the distribution
unique_probes <- names(which(table(affyids$Accession) == 1))
affyids_unique <- affyids[Accession %in% unique_probes]

# select only expression columns relating to prostate cell lines and
# melt data into 2 columns, one for the cell line, one for the expression value
prostate_idx <- grep("PROSTATE", colnames(ccle))
prostate_idx <- c(1, 2, prostate_idx)

# assign HUGO gene names to probe IDs and expression
ccle_annotated <- merge(
    x = ccle[, ..prostate_idx],
    y = affyids_unique,
    all.x = TRUE
)

pca_exprs <- melt(
    ccle_annotated,
    id.vars = c(11, 1, 2)  # melt based on first two columns
)
colnames(pca_exprs) <- c("Gene", "Accession", "Description", "Cell", "Expression")


# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = pca_exprs, aes(x = Cell, y = Expression, fill = Cell))
    + geom_boxplot()
    + geom_point(
        data = pca_exprs[Gene == "FOXA1"],
        mapping = aes(x = Cell, y = Expression),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Cell Line", y = "Gene Expression (RMA-normalized)")
    + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    + guides(fill = FALSE)
)
ggsave(
    "cell-lines-microarray.png",
    height = 12,
    width = 20,
    units = "cm",
    bg = "transparent"
)
