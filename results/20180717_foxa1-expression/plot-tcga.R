# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

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

tcga_percentile <- fread(
    "tcga-expression-percentiles.tsv",
    header = TRUE,
    sep = "\t"
)

# ==============================================================================
# Analysis
# ==============================================================================
# get FOXA1 expression percentiles
foxa1_percentiles <- data.table(
    PatientID = colnames(tcga_percentile[Description == "FOXA1", .SD, .SDcols = 1:550]),
    Percentile = as.vector(t(tcga_percentile[Description == "FOXA1", .SD, .SDcols = 1:550]))
)
# order by descreasing percentile value
foxa1_percentiles <- foxa1_percentiles[order(-rank(Percentile))]
foxa1_percentiles[, PatientID := factor(PatientID, levels = PatientID, ordered = TRUE)]

# melt data together for parsing by ggplot
pca_exprs <- melt(
    tcga,
    id.vars = c(1, 552)  # EnsemblID and HUGO name
)
colnames(pca_exprs) <- c("EnsemblID", "Description", "PatientID", "FPKM")
pca_exprs[, logFPKM := log10(FPKM + 1)]

# melt percentile data
pca_percentile <- melt(
    tcga_percentile,
    id.vars = c(551, 552)  # EnsemblID and HUGO name
)
colnames(pca_percentile) <- c("EnsemblID", "Description", "PatientID", "Percentile")


# ==============================================================================
# Plots
# ==============================================================================
# FPKM plots
gg <- (
    ggplot(data = pca_exprs, aes(x = PatientID, y = logFPKM, fill = PatientID))
    + geom_boxplot()
    + geom_point(
        data = pca_exprs[Description == "FOXA1"],
        mapping = aes(y = logFPKM),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Patient", y = "log10(FPKM + 1)")
    + guides(fill = FALSE)
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1),
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
    "tcga-FPKM.png",
    height = 12,
    width = 100,
    units = "cm",
    bg = "transparent"
)

# FOXA1 percentile plot
gg <- (
    ggplot(
        data = foxa1_percentiles,
        aes(x = PatientID, y = 100 * Percentile)
    )
    + geom_col()
    + labs(x = "Patients", y = "Percentile of FOXA1 Expression")
    + guides(fill = FALSE)
    + coord_cartesian(ylim = c(50, 100))
    + theme(
        # font sizes for axes and legend
        axis.text.x = element_blank(),
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
    "tcga-percentile.png",
    height = 12,
    width = 40,
    units = "cm",
    bg = "transparent"
)
