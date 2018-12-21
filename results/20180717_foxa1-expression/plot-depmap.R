# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))
suppressMessages(library("grid"))
suppressWarnings(library("stringr"))

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
    DepMap_ID
]
#   remove columns corresponding to all non-cancerous lines
noncancer_col_idx = which(depmap[, Line] %in% noncancer_lines)
depmap = depmap[-noncancer_col_idx]

cat("Cleaning metadata\n")
# clean up text for simpler parsing
#   find index of underscore in string
metadata[, Underscore := as.vector(regexpr("_", CCLE_Name))]
#   classify cell line based on tissue type
metadata[, Tissue := substring(CCLE_Name, Underscore + 1, 100)]
metadata[, Tissue := factor(str_to_title(Tissue))]
#   extract readable name of cell line
metadata[, Name := substring(CCLE_Name, 1, Underscore - 1)]
metadata[, Name := factor(Name)]
#   further underscores in Tissue replaced with spaces
metadata[, Tissue := gsub("_", " ", Tissue)]
#   remove unused column
metadata[, Underscore := NULL]

# map cell names to more visually-friendly names, consistent with other plots
metadata[Name == "22RV1", Name := "22Rv1"]
metadata[Name == "LNCAPCLONEFGC", Name := "LNCaP"]
metadata[Name == "MDAPCA2B", Name := "MDA PCa 2B"]
metadata[Name == "NCIH660", Name := "NCI-H660"]
metadata[Name == "PRECLH", Name := "PrEC LH"]
metadata[Name == "VCAP", Name := "VCaP"]

cat("Merging expression data and metadata\n")
# add cell line names to expression data
depmap = merge(
    depmap,
    metadata[, .(DepMap_ID, Name, Tissue)],
    by.x = "Line",
    by.y = "DepMap_ID"
)

# melt into format for plotting
exprs = melt(
    depmap,
    id.vars = c("Line", "Name", "Tissue"),
    variable.name = "Gene",
    value.name = "log2TPM"
)

# select prostate lines
prostate_lines = metadata[
    get("Primary Disease") == "Prostate Cancer",
    .(DepMap_ID, Name)
]
colnames(prostate_lines) = c("ID", "Name")

# get selection of prostate-specific cell lines
pca_exprs = melt(
    depmap[Line %in% prostate_lines[, ID], .SD],
    id.vars = c("Line", "Name", "Tissue"),
    variable.name = "Gene",
    value.name = "log2TPM"
)

# get FOXA1 expression in all lines
foxa1_col_idx = grep("FOXA1", colnames(depmap))
foxa1_exprs = depmap[
    , .SD,
    #           Line, FOXA1, Name, Tissue
    .SDcols = c(1, foxa1_col_idx, dim(depmap)[2] - 1,  dim(depmap)[2])
]
colnames(foxa1_exprs) = c("Line", "FOXA1", "Name", "Tissue")
# get FOXA1 expression in prostate lines
pca_foxa1_exprs = pca_exprs[grep("FOXA1", Gene), .(Name, log2TPM)]


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
        data = pca_foxa1_exprs,
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

# get x-axis ordering for plots
foxa1_order = foxa1_exprs[, median(FOXA1), by = Tissue]
foxa1_exprs[, Tissue := factor(
    Tissue,
    levels = foxa1_order[rev(order(V1)), Tissue],
    ordered = TRUE
)]

gg_count = (
    ggplot(data = foxa1_exprs)
    + geom_bar(aes(x = Tissue), fill = "#FDA328", colour = "black")
    + labs(x = NULL, y = "Number of Cell Lines")
    + theme_classic()
    + theme(
        axis.text.x = element_blank()
    )
)
grob_count = ggplotGrob(gg_count)

gg_exprs = (
    ggplot(data = foxa1_exprs)
    + geom_boxplot(aes(x = Tissue, y = FOXA1), fill = "#EF5350")
    + labs(x = NULL, y = "FOXA1 Expression log2(TPM + 1)")
    + theme_classic()
    + theme(
        axis.text.x = element_text(angle = 60, hjust = 1)
    )
)
grob_exprs = ggplotGrob(gg_exprs)

g = rbind(grob_count, grob_exprs, size = 'first')
g$widths = unit.pmax(grob_count$widths, grob_exprs$widths)

pdf(
    "foxa1-all-lines.pdf",
    height = 12,
    width = 20
)
grid.newpage()
grid.draw(g)
dev.off()
