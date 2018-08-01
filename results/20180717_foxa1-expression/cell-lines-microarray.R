# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

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
# load Affymetrix probe IDs and information
affyids <- fread(
    "Affymetrix_U133_Plus_2-probe-IDs.tsv",
    sep = "\t",
    header = TRUE,
    select = 1:2,
    col.names = c("Accession", "Gene")
)

# don't annotate probe with multiple annotations
# doing so miscounts the number of probes in total, and skews the distributions
unique_probes <- names(which(table(affyids$Accession) == 1))
affyids_unique <- affyids[Accession %in% unique_probes]

# select only expression columns relating to prostate cell lines and
prostate_idx <- grep("PROSTATE", colnames(ccle))
# drop PRECLH_PROSTATE example (can't find any information about it)
prostate_idx <- prostate_idx[-8]
prostate_idx <- c(1, 2, prostate_idx)

# assign HUGO gene names to probe IDs and expression
# only probes with unique HUGO annotations get annotated
# no probes are removed from the actual data
ccle_annotated <- merge(
    x = ccle[, ..prostate_idx],
    y = affyids_unique,
    all.x = TRUE
)

# melt data into 2 columns, one for the cell line, one for the expression value
pca_exprs <- melt(
    ccle_annotated,
    id.vars = c(10, 1, 2)  # melt based on annotations
)
colnames(pca_exprs) <- c("Gene", "Accession", "Description", "Cell", "Expression")

# map cell names to more visually-friendly names
pca_exprs[Cell == "22RV1_PROSTATE", Cell := "22Rv1"]
pca_exprs[Cell == "DU145_PROSTATE", Cell := "DU145"]
pca_exprs[Cell == "LNCAPCLONEFGC_PROSTATE", Cell := "LNCaP"]
pca_exprs[Cell == "MDAPCA2B_PROSTATE", Cell := "MDA PCa 2B"]
pca_exprs[Cell == "NCIH660_PROSTATE", Cell := "NCI-H660"]
pca_exprs[Cell == "PC3_PROSTATE", Cell := "PC3"]
pca_exprs[Cell == "VCAP_PROSTATE", Cell := "VCaP"]

# order the factor for plotting
pca_exprs$Cell <- factor(
    pca_exprs$Cell,
    levels = c(
        "22Rv1",
        "DU145",
        "LNCaP",
        "MDA PCa 2B",
        "NCI-H660",
        "PC3",
        "VCaP"
    ),
    ordered = TRUE
)

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
    + guides(fill = FALSE)
)
ggsave(
    "cell-lines-microarray.png",
    height = 12,
    width = 20,
    units = "cm",
    bg = "transparent"
)
