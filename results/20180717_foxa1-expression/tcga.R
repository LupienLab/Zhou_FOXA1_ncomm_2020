# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read TCGA manifest data
manifest <- fread(
    "../../data/external/TCGA-PRAD/gdc_manifest_20180717_162017.txt",
    header = TRUE,
    sep = "\t"
)

# take first few letters as shorter ID for each patient
# I've checked, and yes these IDs are unique
separate_ids <- strsplit(manifest$id, "-")
small_ids <- sapply(separate_ids, function(x) x[1])
manifest[, PatientID := small_ids]

# read first patient file
tcga <- fread(
    paste0("zcat ../../data/external/TCGA-PRAD/", manifest[1, filename]),
    header = FALSE,
    col.names = c("EnsemblID", manifest[1, PatientID]),
    sep = "\t"
)

# read and merge subsequent patients' data
for (i in 2:manifest[, .N - 1]) {
    print(i)
    newdf = fread(
        paste0("zcat ../../data/external/TCGA-PRAD/", manifest[i, filename]),
        header = FALSE,
        col.names = c("EnsemblID", manifest[i, PatientID]),
        sep = "\t"
    )
    tcga <- merge(tcga, newdf)
}


# ==============================================================================
# Analysis
# ==============================================================================
# strip EnsemblIDs of periods (i.e. their version number, just keep the gene ID)
# all EnsemblIDs are 15 characters long
# see https://useast.ensembl.org/Help/Faq?id=488
tcga[, EnsemblID := substr(EnsemblID, 1, 15)]

# load Ensembl gene IDs
ensembl <- fread(
    "ensembl-gene-IDs.tsv",
    header = TRUE,
    sep = "\t",
    select = 1:2
)

# map Ensembl ID to a gene name
tcga <- merge(
    x = tcga,
    y = ensembl,
    by.x = "EnsemblID",
    by.y = "ensembl_gene_id",
    all.x = TRUE
)

# find FOXA1 expression percentile wrt all other mapped reads for that patient
tcga_percentile <- as.data.table(apply(
    tcga[, 2:11],
    2,
    function(x) trunc(rank(x, na.last = NA))/sum(!is.na(x))
))
tcga_percentile[, EnsemblID := tcga$EnsemblID]
tcga_percentile[, Description := tcga$hgnc_symbol]
tcga_percentile[Description == "FOXA1"]

# melt data together for parsing by ggplot
pca_exprs <- melt(
    tcga,
    id.vars = c(1, 552)  # EnsemblID and HUGO name
)
# pca_exprs <- melt(
#     tcga2[, c(1:100, 552)],
#     id.vars = c(1, 101)  # EnsemblID and HUGO name
# )
colnames(pca_exprs) <- c("EnsemblID", "Description", "PatientID", "FPKM")
pca_exprs[, logFPKM := log10(FPKM + 1)]



# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = pca_exprs, aes(x = PatientID, y = logFPKM, fill = PatientID))
    + geom_boxplot()
    + geom_point(
        data = pca_exprs[Description == "FOXA1"],
        mapping = aes(y = logFPKM),
        fill = "red",
        pch = 21,
        size = 4
    )
    + labs(x = "Patient", y = "log10(FPKM + 1)")
    + guides(fill = FALSE)
    + theme(
        # font sizes for axes and legend
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1)
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
    "tcga.png",
    height = 12,
    width = 100,
    units = "cm",
    bg = "transparent"
)
