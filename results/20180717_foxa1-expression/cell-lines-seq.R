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
prostate_idx <- grep("PROSTATE", colnames(depmap))
# drop PRECLH_PROSTATE example (can't find any information about it)
prostate_idx <- prostate_idx[-7]
prostate_idx <- c(1, 2, prostate_idx)
# melt data into 2 columns, one for the cell line, one for the expression value
pca_exprs <- melt(
    depmap[, ..prostate_idx],  # reduce to only columns of interest
    id.vars = c(1, 2)  # melt based on first two columns
)
colnames(pca_exprs) <- c("Name", "Description", "Cell", "Expression")

# map cell names to more visually-friendly names
pca_exprs[Cell == "22RV1_PROSTATE", Cell := "22Rv1"]
pca_exprs[Cell == "DU145_PROSTATE", Cell := "DU145"]
pca_exprs[Cell == "LNCAPCLONEFGC_PROSTATE", Cell := "LNCaP"]
pca_exprs[Cell == "MDAPCA2B_PROSTATE", Cell := "MDA PCa 2B"]
pca_exprs[Cell == "NCIH660_PROSTATE", Cell := "NCI-H660"]
pca_exprs[Cell == "PC3_PROSTATE", Cell := "PC3"]
pca_exprs[Cell == "VCAP_PROSTATE", Cell := "VCaP"]

# order the factor
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
        mapping = aes(x = Cell, y = log2(Expression + 1)),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Cell Line", y = "Gene Expression log2(RPKM + 1)")
    + guides(fill = FALSE)
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        # plot background colouring
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#9e9e9e"),
        panel.background = element_rect(fill = "transparent")
    )
)
ggsave(
    "cell-lines-seq.png",
    height = 12,
    width = 20,
    units = "cm",
    bg = "transparent"
)
