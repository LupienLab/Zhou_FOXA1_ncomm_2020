# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("ggplot2"))
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read in preprocessed data
all_table <- fread(
    "essentiality-scores.tsv",
    header = TRUE,
    sep = "\t"
)

crispr <- all_table[Method == "CRISPR"]
# crispr[, Tissue := factor(Tissue)]
rnai <- all_table[Method == "RNAi"]
# rnai[, Tissue := factor(Tissue)]

# ==============================================================================
# Plots
# ==============================================================================
gg_essentiality <- (
    ggplot(
        data = all_table[Tissue == "Prostate"],
        mapping = aes(x = Line, y = Score, fill = Method)
    )
    + geom_boxplot()
    + geom_point(
        data = all_table[Description == "FOXA1" & Tissue == "Prostate"],
        mapping = aes(x = Line, y = Score),
        fill = "red",
        pch = 21,
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
    filename = "essentiality-prostate-all.png",
    plot = gg_essentiality,
    height = 12,
    width = 20,
    units = "cm"
)

# CRISPR knockout plots for FOXA1 in all cell lines
gg_foxa1_crispr <- (
    ggplot(data = crispr)
    + geom_boxplot(
        aes(x = Line, y = Score, fill = Tissue)
    )
    + geom_point(
        data = crispr[Description == "FOXA1"],
        mapping = aes(x = Line, y = Score),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = FALSE)
    + ggtitle("CRISPR Essentiality Scores")
    + facet_grid(. ~ Tissue, scales = "free_x")
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 6, angle = 90, vjust = 1),
        strip.text.x = element_text(size = 6),
        # plot background colouring
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#9e9e9e"),
        panel.background = element_rect(fill = "transparent")
    )
)
ggsave(
    filename = "essentiality-all-crispr.png",
    plot = gg_foxa1_crispr,
    height = 12,
    width = 40,
    units = "cm"
)

# RNAi plots for FOXA1 in all cell lines
gg_foxa1_rnai <- (
    ggplot(data = rnai)
    + geom_boxplot(
        aes(x = Line, y = Score, fill = Tissue)
    )
    + geom_point(
        data = rnai[Description == "FOXA1"],
        mapping = aes(x = Line, y = Score),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = FALSE)
    + ggtitle("RNAi Essentiality Scores")
    + facet_wrap(~ Tissue, scales = "free_x")
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 6, angle = 90, vjust = 1),
        strip.text.x = element_text(size = 8),
        # plot background colouring
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#9e9e9e"),
        panel.background = element_rect(fill = "transparent")
    )
)
ggsave(
    filename = "essentiality-all-rnai.png",
    plot = gg_foxa1_rnai,
    height = 30,
    width = 30,
    units = "cm"
)
