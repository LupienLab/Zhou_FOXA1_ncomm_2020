# =================================================================================================
# Environment
# =================================================================================================
library("ggplot2")
library("data.table")

# =================================================================================================
# Data
# =================================================================================================
cell_proliferation <- fread(
    "../../data/cell-proliferation.tsv",
    header = TRUE
)

# =================================================================================================
# Plot
# =================================================================================================
gg_proliferation <- (
    # data
    ggplot(
        data = cell_proliferation,
        aes(
            x = Day,
            y = Count / 1e3,
            fill = Method,
            colour = Method,
            group = Method
        )
    )
    # dots
    + geom_point(
        colour = "black",
        shape = 21,
        size = 4,
        stroke = 1
    )
    # lines
    + stat_summary(
        fun.data = mean_se,
        geom = "line",
        size = 2,
        alpha = 0.75
    )
    # error bars
    + stat_summary(
        fun.data = mean_se,
        geom = "errorbar",
        width = 0.5,
        size = 1,
        alpha = 0.75
    )
    + labs(y = "Cell Count (1000s)")
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        # plot background colouring
        panel.grid.major = element_line(colour = "#757575"),
        panel.grid.minor = element_line(colour = "#AAAAAA"),
        panel.background = element_rect(fill = "transparent")
    )
)

# save full plot
ggsave(
    gg_proliferation,
    file = "cell-proliferation_0-5.png"
)

# save plot between starting at day 2
ggsave(
    gg_proliferation + lims(x = c(1.75, 5.25)),
    file = "cell-proliferation_2-5.png"
)