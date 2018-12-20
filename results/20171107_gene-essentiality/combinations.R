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

prostate_lines = rnai[Tissue == "Prostate", unique(Line)]
foxa1 = rnai[Gene == "FOXA1"]
all_lines = foxa1[!is.na(Score), unique(Line)]

# ==============================================================================
# Analysis
# ==============================================================================
# Test: For m randomly selected lines, what is the probability that
# median({e_i}) < median({e_prostate})
prostate_median = foxa1[Line %in% prostate_lines, median(Score)]

# number of combinations for each method
n_perm_possible = choose(foxa1[, .N], length(prostate_lines))
n_perm = 1e6

# RNAi combinations
cat("Calculating permutations\n")
medians = rep(NA, n_perm)
for (i in 1:n_perm) {
    if (i %% 1e3 == 0) cat(i, "\n")
    sampled_lines = sample(all_lines, 8, replace = FALSE)
    medians[i] = foxa1[Line %in% sampled_lines, median(Score)]
}

p = (length(which(medians == prostate_median)) + 1) / n_perm

p_vals = data.table(
    N = n_perm,
    P = p
)

fwrite(
    p_vals,
    "essentiality-combinations.tsv",
    col.names = TRUE,
    sep = "\t"
)
