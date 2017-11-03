# =================================================================================================
# Environment
# =================================================================================================
library("ggplot2")
library("data.table")

# =================================================================================================
# Data
# =================================================================================================
mutated_regions <- fread(
    "../../data/processed/Mutated_Regions.txt",
    header = TRUE
)
# don't consider DistalREs
# mutated_promoters <- mutated_regions[peakANNO != "DistalRE"]
# mutated_promoters <- mutated_regions
# reorder for bar plot
# mutated_promoters$peakID <- reorder(mutated_promoters$peakID, -mutated_promoters$neglog10qval)

# significance level
sig_level <- 0.01

# select top candidates to include in inset plot
n_sig <- length(which(mutated_promoters[, neglog10qval >= -log10(sig_level)]))

# =================================================================================================
# Plots
# =================================================================================================
# full screen plot
gg_mutated_full <- (
    ggplot(
        data = mutated_regions,
        aes(
            x = peakMUTr,
            y = neglog10qval
        )
    )
    + geom_point()
    + labs(
        x = "Mutation Rate",
        y = "-log10(FDR)"
    )
    + scale_x_discrete(breaks = NULL)
)
print(gg_mutated_full)
# save full plot
ggsave(
    gg_mutated_full,
    file = "smurf_full.png"
)


# top regions, inset plot
gg_mutated_top <- (
    # data of most significantly mutated regions
    ggplot(
        data = mutated_promoters[1:n_sig]
    )
    # column plot
    + geom_col(
        aes(
            x = peakID,
            y = neglog10qval
        )
    )
    # change labels on x-axis to name of associated gene
    # + scale_x_discrete(
    #     labels = mutated_promoters[1:66]$peakANNO
    # )
    # axes labels
    + labs(
        x = "Gene-associated Promoter",
        y = "-log10(FDR)"
    )
    + theme(
        # rotate and adjust x-axis names
        axis.text.x = element_text(angle = 90, hjust = 1),
        # font sizes for axes and legend
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 16),
        # plot background colouring
        axis.ticks = element_line(colour = NA)
    )
)

# save full plot
ggsave(
    gg_mutated_top,
    file = "smurf_top.png"
)