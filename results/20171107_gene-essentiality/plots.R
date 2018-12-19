# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("ggplot2"))
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read in preprocessed data
rnai = fread(
    "essentiality-scores.tsv",
    header = TRUE,
    sep = "\t"
)


foxa1 = rnai[Gene == "FOXA1"]
foxa1_medians = foxa1[, median(Score, na.rm = TRUE), by = Tissue]
foxa1[, Tissue := factor(Tissue, levels = foxa1_medians[order(V1), Tissue], ordered = TRUE)]

# ==============================================================================
# Plots
# ==============================================================================
gg_essentiality = (
    ggplot(
        data = rnai[Tissue == "Prostate"],
        mapping = aes(x = Line, y = Score)
    )
    + geom_violin(fill = "#26A69A")
    + geom_boxplot(width = 0.05, fill = "#26A69A")
    + geom_point(
        data = rnai[Gene == "FOXA1" & Tissue == "Prostate"],
        mapping = aes(x = Line, y = Score),
        fill = "#EF5350",
        pch = 21,
        size = 4,
        alpha = 0.8
    )
    + geom_point(
        data = rnai[Gene == "FOXA1" & Tissue == "Prostate"],
        mapping = aes(x = Line, y = Score),
        fill = "black",
        pch = 21,
        size = 0.25
    )
    + labs(x = NULL, y = "Gene Essentiality Score")
    + guides(fill = FALSE)
    + scale_y_continuous(breaks = seq(-3, 3, 0.5), limits = c(-2.7, 2))
    + theme_classic()
    + theme(
        # font sizes for axes and legend
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0),
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
    filename = "essentiality-prostate.png",
    plot = gg_essentiality,
    height = 12,
    width = 20,
    units = "cm"
)

# RNAi plots for FOXA1 in all cell lines compared to all other genes
gg_all = (
    ggplot()
    + geom_boxplot(
        data = rnai,
        mapping = aes(x = Line, y = Score, fill = Tissue)
    )
    + geom_point(
        data = foxa1,
        mapping = aes(x = Line, y = Score),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = FALSE)
    + facet_wrap(~ Tissue, scales = 'free_x')
    + theme_classic()
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 8, angle = 90, vjust = 1),
        strip.text.x = element_text(size = 8),
        # plot background colouring
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#9e9e9e"),
        panel.background = element_rect(fill = "transparent")
    )
)
ggsave(
    filename = "essentiality-all.png",
    plot = gg_all,
    height = 30,
    width = 30,
    units = "cm"
)

# RNAi plots for FOXA1 in all cell lines, grouped by tissue
gg_foxa1 = (
    ggplot(
        data = foxa1,
        mapping = aes(x = Tissue, y = Score)
    )
    + stat_boxplot(geom ='errorbar')
    + geom_boxplot(fill = "#26A69A")
    + geom_point(
        position = position_jitter(width = 0.1)
    )
    + labs(x = NULL, y = "FOXA1 Essentiality Score")
    + theme_classic()
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 6, angle = 90, hjust = 1),
        strip.text.x = element_text(size = 8),
        # plot background colouring
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#9e9e9e"),
        panel.background = element_rect(fill = "transparent")
    )
)
ggsave(
    filename = "essentiality-foxa1.png",
    plot = gg_foxa1,
    height = 20,
    width = 20,
    units = "cm"
)
