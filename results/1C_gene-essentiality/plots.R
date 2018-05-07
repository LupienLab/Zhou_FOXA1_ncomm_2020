# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("ggplot2"))
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# Achilles data
achilles_crispr <- fread(
    "../../data/external/Achilles_v3.3.8.Gs.gct",
    header = TRUE,
    skip = 2
)
achilles_rnai <- fread(
    "../../data/external/Achilles_QC_v2.4.3.rnai.Gs.gct",
    header = TRUE,
    skip = 2
)

# Gene of interest
gene <- "FOXA1"

# Cell lines of interest
crispr_lines <- grep("PROSTATE", names(achilles_crispr))
rnai_lines <- grep("PROSTATE", names(achilles_rnai))

print(names(achilles_crispr)[crispr_lines])
# > [1] "LNCAPCLONEFGC_PROSTATE" "PC3_PROSTATE"

print(names(achilles_rnai)[rnai_lines])
# > [1] "22RV1_PROSTATE"   "NCIH660_PROSTATE" "VCAP_PROSTATE"

# transforming data for plotting
achilles_all <- data.table(
    Description = c(
        rep(achilles_crispr[, Description], length(crispr_lines)),
        rep(achilles_rnai[, Description], length(rnai_lines))
    ),
    Dataset = factor(c(
        rep("CRISPR", nrow(achilles_crispr) * length(crispr_lines)),
        rep("RNAi", nrow(achilles_rnai) * length(rnai_lines))
    )),
    Cell = factor(c(
        rep("LNCaP", nrow(achilles_crispr)),
        rep("PC3", nrow(achilles_crispr)),
        rep("22RV1", nrow(achilles_rnai)),
        rep("NCIH660", nrow(achilles_rnai)),
        rep("VCaP", nrow(achilles_rnai))
    )),
    Score = c(
        achilles_crispr[, LNCAPCLONEFGC_PROSTATE],
        achilles_crispr[, PC3_PROSTATE],
        achilles_rnai[, `22RV1_PROSTATE`],
        achilles_rnai[, NCIH660_PROSTATE],
        achilles_rnai[, VCAP_PROSTATE]
    )
)

# melted data.table for essentiality across all cell lines
crispr_table <- melt(
    achilles_crispr,
    id.vars = "Description",
    measure.vars = 3:35,
    variable.name = "Cell",
    value.name = "Score"
)
# find index of underscore in string
crispr_table[, Underscore := as.vector(regexpr("_", Cell))]
# classify cell line based on tissue type
crispr_table[, Type := substring(Cell, Underscore + 1, 100)]

# melted data.table for essentiality across all cell lines
rnai_table <- melt(
    achilles_rnai,
    id.vars = "Description",
    measure.vars = 3:218,
    variable.name = "Cell",
    value.name = "Score"
)
# find index of underscore in string
rnai_table[, Underscore := as.vector(regexpr("_", Cell))]
# classify cell line based on tissue type
rnai_table[, Type := substring(Cell, Underscore + 1, 100)]

# ==============================================================================
# Plots
# ==============================================================================
gg_essentiality <- (
    ggplot(
        data = achilles_all,
        mapping = aes(x = Cell, y = Score, fill = Dataset)
    )
    + geom_boxplot()
    + geom_point(
        data = achilles_all[Description == "FOXA1"],
        mapping = aes(x = Cell, y = Score),
        fill = "red",
        pch = 23,
        size = 4
    )
    + labs(x = "Cell Line", y = "Gene Essentiality Score")
    + guides(fill = guide_legend(title = "Method"))
    + lims(y = c(-5, 5))
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
    gg_essentiality,
    file = "essentiality.png"
)

# CRISPR knockout plots for FOXA1 in all cell lines
gg_foxa1_crispr <- (
    ggplot(data = crispr_table)
    + geom_boxplot(
        aes(x = Cell, y = Score, fill = Type)
    )
    + geom_point(
        data = crispr_table[Description == "FOXA1"],
        mapping = aes(x = Cell, y = Score),
        fill = "red",
        pch = 23,
        size = 4
    )
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = FALSE)
    + ggtitle("CRISPR Essentiality Scores")
    + facet_grid(. ~ Type, scales = "free_x")
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 6, angle = 90, vjust = 1),
    )
)
ggsave(
    "crispr-essentiality-all.png",
    height = 18,
    width = 40,
    units = "cm"
)

# RNAi plots for FOXA1 in all cell lines
gg_foxa1_rnai <- (
    ggplot(data = rnai_table)
    + geom_boxplot(
        aes(x = Cell, y = Score, fill = Type)
    )
    + geom_point(
        data = rnai_table[Description == "FOXA1"],
        mapping = aes(x = Cell, y = Score),
        fill = "red",
        pch = 23,
        size = 4
    )
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = FALSE)
    + ggtitle("RNAi Essentiality Scores")
    + facet_grid(. ~ Type, scales = "free_x")
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 6, angle = 90, vjust = 1),
    )
)
ggsave(
    "rnai-essentiality-all.png",
    height = 18,
    width = 40,
    units = "cm"
)
