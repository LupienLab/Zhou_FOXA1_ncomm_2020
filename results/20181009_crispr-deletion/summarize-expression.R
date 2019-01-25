# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# load CRISPR data
crispr_rel_noguides = fread(
    "CRISPR_NormToNoGuides.all.tsv",
    sep = "\t",
    fill = TRUE,
    col.names = c("Sample", "Expression", "Replicate")
)

# load lentiviral data
lentiviral_rel_chr14 = fread(
    "Lentiviral_Cas9.all.tsv",
    sep = "\t",
    fill = TRUE,
    col.names = c("Sample", "Expression", "Replicate")
)

# ==============================================================================
# Analysis
# ==============================================================================
summary_crispr = crispr_rel_noguides[, .(mean(Expression), sd(Expression)), by = Sample]
colnames(summary_crispr) = c("Sample", "Expression", "SD")

summary_lentiviral = lentiviral_rel_chr14[, .(mean(Expression), sd(Expression)), by = Sample]
colnames(summary_lentiviral) = c("Sample", "Expression", "SD")

# ==============================================================================
# Output
# ==============================================================================
# save summary tables
fwrite(
    summary_crispr,
    "CRISPR_NormToNoGuides.summary.tsv",
    sep = '\t',
    col.names = TRUE
)
fwrite(
    summary_lentiviral,
    "Lentiviral_Cas9.summary.tsv",
    sep = '\t',
    col.names = TRUE
)
