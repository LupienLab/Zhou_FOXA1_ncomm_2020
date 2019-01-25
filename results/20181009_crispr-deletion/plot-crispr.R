# ==============================================================================
# Environment
# ==============================================================================
library("ggplot2")
library("GenomicRanges")
library("data.table")
library("ggpubr")
library("ggsignif")
library("dplyr")

# ==============================================================================
# Data
# ==============================================================================
expn_rel_noguides = fread(
    "CRISPR_NormToNoGuides.all.tsv",
    sep = "\t",
    fill = TRUE,
    header = TRUE,
    col.names = c("Sample", "Expression", "Replicate")
)

summary_rel_noguides = fread(
    "CRISPR_NormToNoGuides.summary.tsv",
    sep = "\t",
    col.names = c("Sample", "Mean", "SD")
)


#==============================================================================
# t-test for hotspot regions regulating FOXA1 expression by CRISPR/Cas9 deletion
# (contrasting to Neg #1)

wilcox.test(
    expn_rel_noguides[Sample == "CRE1", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "CRE2", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "CRE4", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "CRE1 + CRE2", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "CRE1 + CRE4", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "CRE2 + CRE4", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "FOXA1 Promoter (+)", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)
wilcox.test(
    expn_rel_noguides[Sample == "AAVS1 (-)", Expression],
    expn_rel_noguides[Sample == "Chr14 (-)", Expression]
)

#==============================================================================
#Plotting data with ggplot2
#specify data
# add "path" object to show dots
#specify x/y axes and labels
#point graph overlapped with bargraph

LNCaP <- (
    ggplot()
    #bargraph of mean
    + geom_bar(
        data = summary_rel_noguides,
        mapping = aes(
            x = factor(
                Sample,
                levels = c(
                    "CRE1", 
                    "CRE2", 
                    "CRE4", 
                    "CRE1 + CRE2", 
                    "CRE1 + CRE4", 
                    "CRE2 + CRE4", 
                    "FOXA1 Promoter (+)", 
                    "AAVS1 (-)", 
                    "Chr14 (-)"
                ),
            ordered = TRUE
            ),
            y = Mean
        ),
        stat = "identity",
        fill = c(
            rep("plum1", 3), 
            rep("darksalmon", 3), 
            rep("cornflowerblue", 1), 
            rep("grey", 2)
        ),
        width = .5,
        position="dodge",
        size = 0.7,
        colour = "black"
    )
    #errorbars of data (standard deviation of mean)
    + geom_errorbar(
        data = summary_rel_noguides,
        mapping = aes(x = Sample, ymin = Mean - SD, ymax = Mean + SD),
        width = 0.2
    )
    #points to showcase datapoints
    + geom_point(
        data = expn_rel_noguides,
        mapping = aes(x = Sample, y = Expression), 
        size = 4, 
        color = "black",
        shape = 18
    )
    + ggtitle("FOXA1 Expression Upon Hotspot Deletion via CRISPR/Cas9")
    + labs(x = NULL, y = "Relative FOXA1 Expression (%)")
    + scale_y_continuous(expand = c(0,0), limit = c(0, 130))
    + theme_bw()
    + theme(
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
    )
)


#==============================================================================
#Modifying font
LNCaP_no_font <- LNCaP + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_line()
)

#==============================================================================
#Saving the figure
ggsave("crispr.png", plot = LNCaP, width = 13, height = 8, dpi = 700)
ggsave("crispr-no-font.png", plot = LNCaP_no_font, width = 13, height = 8, dpi = 700)
