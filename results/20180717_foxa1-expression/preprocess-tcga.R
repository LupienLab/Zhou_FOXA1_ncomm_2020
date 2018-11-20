# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

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
for (i in 2:manifest[, .N]) {
    cat(i, "of", manifest[, .N], "\n")
    newdf = fread(
        paste0("zcat ../../data/external/TCGA-PRAD/", manifest[i, filename]),
        header = FALSE,
        col.names = c("EnsemblID", manifest[i, PatientID]),
        sep = "\t"
    )
    tcga <- merge(tcga, newdf)
}


# ==============================================================================
# Save output
# ==============================================================================
# strip EnsemblIDs of periods (i.e. their version number, just keep the gene ID)
# all EnsemblIDs are 15 characters long
# see https://useast.ensembl.org/Help/Faq?id=488
tcga[, EnsemblID := substr(EnsemblID, 1, 15)]
# save data for future use
fwrite(
    tcga,
    "tcga-aggregated.tsv",
    col.names = TRUE,
    sep = "\t"
)
