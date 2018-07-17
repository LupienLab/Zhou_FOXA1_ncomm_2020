# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read depmap data
depmap <- fread(
    "../../data/external/CCLE/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct",
    header = TRUE,
    sep = "\t",
    skip = 2  # skip "#1.2" and dimension lines
)

# select only expression columns relating to prostate cell lines and
# melt data into 2 columns, one for the cell line, one for the expression value
prostate_idx <- grep("PROSTATE", colnames(depmap))
prostate_idx <- c(1, 2, prostate_idx)
pca_exprs <- melt(
    depmap[, ..prostate_idx],  # reduce to only columns of interest
    id.vars = c(1, 2)  # melt based on first two columns
)
colnames(pca_exprs) <- c("Name", "Description", "Cell", "Expression")

# ==============================================================================
# Analysis
# ==============================================================================
# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = pca_exprs, aes(x = Cell, y = log10(Expression + 1), fill = Cell))
    + geom_violin()
    + geom_point(
        data = pca_exprs[Description == "FOXA1"],
        mapping = aes(x = Cell, y = log10(Expression + 1)),
        fill = "red",
        pch = 23,
        size = 4
    )
    + labs(x = "Cell Line", y = "log10(RPKM + 1)")
)
ggsave(
    "cell-lines-seq.png",
    height = 25,
    width = 40,
    units = "cm",
    bg = "transparent"
)
