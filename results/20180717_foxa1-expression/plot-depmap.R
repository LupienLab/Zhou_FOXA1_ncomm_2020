# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
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
foxa1_exprs = pca_exprs[grep("FOXA1", Gene), .(Name, log2TPM)]


# ==============================================================================
# Plots
# ==============================================================================
gg = (
    ggplot(
        data = pca_exprs,
        mapping = aes(x = Name, y = log2TPM)
    )
    + geom_boxplot(fill = "#26A69A")
    + geom_point(
        data = pca_exprs[grep("FOXA1", Gene), .SD],
        mapping = aes(x = Name, y = log2TPM),
        fill = "#EF5350",
        pch = 21,
        size = 4
    )
    + labs(x = NULL, y = "Gene Expression log2(TPM + 1)")
    + guides(fill = FALSE)
    + theme_classic()
    + theme(
        # font sizes for axes and legend
        axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 12),
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
    "cell-lines.pdf",
    height = 12,
    width = 20,
    units = "cm",
    bg = "transparent"
)
