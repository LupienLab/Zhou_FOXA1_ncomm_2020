# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read TCGA data
tcga = readRDS("../../data/external/TCGA-PRAD/TCGA_Xena_PRAD.rds")
tcga_dt = as.data.table(tcga)
tcga_dt[, Gene := rownames(tcga)]

percentiles <- fread(
    "tcga-expression-percentiles.tsv",
    header = TRUE,
    sep = "\t"
)
colnames(percentiles)[1] = "Gene"

# ==============================================================================
# Analysis
# ==============================================================================
# get FOXA1 expression percentiles
foxa1_percentiles <- data.table(
    PatientID = colnames(percentiles[Gene == "FOXA1", .SD, .SDcols = 2:498]),
    Percentile = as.vector(t(percentiles[Gene == "FOXA1", .SD, .SDcols = 2:498]))
)
# order by descreasing percentile value
foxa1_percentiles <- foxa1_percentiles[order(-rank(Percentile))]
foxa1_percentiles[, PatientID := factor(PatientID, levels = PatientID, ordered = TRUE)]
# print expression percentiles
print(foxa1_percentiles[Percentile > 0.95, .N])
print(foxa1_percentiles[Percentile > 0.95, .N] / foxa1_percentiles[, .N])

# melt data together for parsing by ggplot
pca_exprs <- melt(
    tcga_dt,
    id.vars = "Gene",
    variable.name = "PatientID",
    value.name = "RSEM"
)

# ==============================================================================
# Plots
# ==============================================================================
# FOXA1 percentile plot
cat("Plotting TCGA percentiles\n")
gg <- (
    ggplot(
        data = foxa1_percentiles,
        aes(x = PatientID, y = 100 * Percentile)
    )
    + geom_col()
    + labs(x = "Patients", y = "Percentile of FOXA1 Expression")
    + guides(fill = FALSE)
    + coord_cartesian(ylim = c(80, 100))
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
    "tcga-percentile-col.png",
    height = 12,
    width = 40,
    units = "cm",
    bg = "transparent"
)

gg <- (
    ggplot(
        data = foxa1_percentiles,
        aes(x = PatientID, y = 100 * Percentile)
    )
    + geom_point()
    + labs(x = "Patients", y = "Percentile of FOXA1 Expression")
    + guides(fill = FALSE)
    + coord_cartesian(ylim = c(80, 100))
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
    "tcga-percentile-point.png",
    height = 12,
    width = 40,
    units = "cm",
    bg = "transparent"
)

# RSEM plots
cat("Plotting RSEM data\n")
gg <- (
    ggplot(data = pca_exprs, aes(x = PatientID, y = RSEM, fill = PatientID))
    + geom_boxplot()
    + geom_point(
        data = pca_exprs[Gene == "FOXA1"],
        mapping = aes(y = RSEM),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Patient", y = "RSEM")
    + guides(fill = FALSE)
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
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
    "tcga-RSEM.png",
    height = 12,
    width = 100,
    units = "cm",
    bg = "transparent"
)
