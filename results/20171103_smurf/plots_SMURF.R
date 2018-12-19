# =================================================================================================
# Environment
# =================================================================================================
library("ggplot2")
library("data.table")

# =================================================================================================
# Data
# =================================================================================================
mutated_regions <- fread(
    "../smurf/Mutated_Regions.txt",
    header = TRUE
)

# significance level
sig_level <- 0.01

# select top candidates to include in inset plot
n_sig <- length(which(mutated_regions[, neglog10qval >= -log10(sig_level)]))

# =================================================================================================
# Plots
# =================================================================================================
# full screen plot
gg_mutated_full <- (
    ggplot(
        data = mutated_regions,
        aes(
            x = peakMUTr,
            y = neglog10qval,
            colour = as.factor(ifelse(peakANNO == "DistalRE", "DistalRE", "Promoter"))
        )
    )
    + geom_point(alpha = 0.5)
    + labs(
        x = "Mutation Rate",
        y = "-log10(FDR)"
    )
    + guides(
        colour = guide_legend(title = "Region Type")
    )
    + geom_hline(
        yintercept = -log10(sig_level),
        linetype = 2
    )
)

# print(gg_mutated_full)
# save full plot
ggsave(
    gg_mutated_full,
    file = "smurf_full.png",
    height = 30,
    width = 40,
    units = "cm"
)

# reorder data for plot
mutated_regions_top <- mutated_regions[1:n_sig]
mutated_regions_top$peakID <- reorder(mutated_regions_top$peakID, - mutated_regions_top$neglog10qval)

# top regions, inset plot
gg_mutated_top <- (
    # data of most significantly mutated regions
    ggplot(
        data = mutated_regions_top
    )
    # column plot
    + geom_col(
        aes(
            x = peakID,
            y = neglog10qval
        )
    )
    # axes labels
    + labs(
        x = "Gene-associated Promoter",
        y = "-log10(FDR)"
    )
    + theme(
        # rotate and adjust x-axis names
        axis.text.x = element_text(angle = 90, hjust = 1),
        # font sizes for axes and legend
        axis.text = element_text(size = 4),
        axis.title = element_text(size = 16),
        # plot background colouring
        axis.ticks = element_line(colour = NA)
    )
)

# print(gg_mutated_top)
# save full plot
ggsave(
    gg_mutated_top,
    file = "smurf_top.png",
    height = 30,
    width = 40,
    units = "cm"
)