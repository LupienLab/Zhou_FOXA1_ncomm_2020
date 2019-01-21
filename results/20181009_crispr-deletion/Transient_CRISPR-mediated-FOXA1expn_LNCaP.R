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
# expression data
expn_rel_noguides <- fread(
    "FOXA1_LNCAP_Cas9_NormToNoGuideControl.tsv",
    sep = "\t",
    fill = TRUE,
    col.names = c("Sample", "Expression", "Replicate")
)

# summary expression data
summary_rel_noguides <- fread(
    "Summary_NormToNoGuidesdata_LNCaP.tsv",
    sep = "\t",
    col.names = c("Sample", "Mean", "SD")
)
#==============================================================================
# T-test for hotspot regions regulating FOXA1 expression by CRISPR/Cas9 deletion
#   (contrasting to Neg#1)
# a is A/e1
# b is B/e2
# c is C/e3
# d is e1+e2
# e is e1+e3
# f is e2+e3
# g is FOXA1 promoter
# h is chr14-
# i is aavs1-


a <- expn_rel_noguides[1:5, 2]
b <- expn_rel_noguides[6:10, 2]
c <- expn_rel_noguides[11:15, 2]
d <- expn_rel_noguides[16:20, 2]
e <- expn_rel_noguides[21:25, 2]
f <- expn_rel_noguides[26:30, 2]
g <- expn_rel_noguides[31:35, 2]
h <- expn_rel_noguides[36:40, 2]
i <- expn_rel_noguides[41:45, 2]


t.test(a,i)
t.test(b,i)
t.test(c,i)
t.test(d,i)
t.test(e,i)
t.test(f,i)
t.test(g,i)
t.test(i,i)

#==============================================================================
#Setting up data for plot

expn_df <- data.table(
#Sample = 
expn_rel_noguides[1:45, 1],
#Expression = 
expn_rel_noguides[1:45, 2])


summary_df <- data.table(
    Sample = expn_rel_noguides[, unique(Sample)],
    Mean = expn_rel_noguides[, mean(Expression), by = Sample]$V1,
    SD = expn_rel_noguides[, sd(Expression), by = Sample]$V1
)

# my_comparisons <- list(
#     c(81.60167, 78.38921, 89.79166, 80.70585, 82.24661),
#     c(100, 100, 100, 100, 100)
# )

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
        aes(
            x = factor(
                Sample,
                levels = c(
                    "CRE1", 
                    "CRE2", 
                    "CRE3", 
                    "CRE1 + CRE2", 
                    "CRE1 + CRE3", 
                    "CRE2 + CRE3", 
                    "FOXA1 Promoter (+)", 
                    "Chr14 (-)", 
                    "AAVS1 (-)"
                ),
                ordered = TRUE
            ),
            y = Mean
        ),
        stat="identity",
        fill = c(
            rep("chocolate1", 1), 
            rep("orchid3", 1),
            rep("seagreen4", 1),
            rep("dodgerblue3", 3), 
            rep("mediumaquamarine", 1), 
            rep("grey", 2)
        ),
        width = .5,
        position = "dodge",
        size = 0.7,
        colour = "black"
    )
    #errorbars of data (standard deviation of mean)
    + geom_errorbar(
        data = summary_df,
        aes(x = Sample, ymin = Mean - SD, ymax = Mean + SD),
        width = 0.2
    )
    #points to showcase datapoints
    + geom_point(
        data = expn_rel_noguides,
        aes(x = Sample, y = Expression), 
        size = 4, 
        color = "black",
        shape = 18
    )
    # + stat_compare_means(
    #   comparisons = my_comparisons
    # )
    # + # Add pairwise comparisons p-value
    #   stat_compare_means(
    #     label.y = 50
    #   )     # Add global p-value
    + ggtitle("FOXA1 Expression Upon Hotspot Deletion via CRISPR/Cas9")
    + labs(y = "Relative FOXA1 Expression (%)")
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
LNCaP_plus_font <- LNCaP + theme(
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 30),
    axis.text.x = element_text(color = "black", size = 24, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 24),
    axis.line = element_line(colour = 'black', size = 1)
)

#==============================================================================
#Saving the figure
ggsave(
    "transient_single_double_dels.png",
    width = 14,
    height = 10,
    dpi = 700
)