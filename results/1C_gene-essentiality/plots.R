# ==============================================================================
# Environment
# ==============================================================================
library("ggplot2")
library("data.table")

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
    + labs(x = "Cell Line", y = "Essentiality Score")
    + guides(fill = guide_legend(title = "Knockdown Method"))
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
        # panel.grid.minor = element_line(colour = "#bdbdbd"),
        panel.background = element_rect(fill = "transparent")
    )
)

ggsave(
    gg_essentiality,
    file = "essentiality.png"
)
