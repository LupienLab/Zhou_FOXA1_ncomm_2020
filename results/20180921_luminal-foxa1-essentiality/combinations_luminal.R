# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("ggplot2"))
suppressMessages(library("data.table"))

# ==============================================================================
# Data
# ==============================================================================
# read in preprocessed data
all_table <- fread(
    "../1C_gene-essentiality/essentiality-scores.tsv",
    header = TRUE,
    sep = "\t"
)

crispr <- all_table[Method == "CRISPR"]
rnai <- all_table[Method == "RNAi"]

# there are 3 cell lines whose classifications all agree between classification methods:
# T-47D, MCF-7, and EFM-19
# T-47D is in the CRISPRi dataset, and the other 2 are in the RNAi dataset

# ==============================================================================
# Analysis
# ==============================================================================
# Test 1: For a randomly selected line, what is the probability that
# e_test < e_luminal = P(e_test < e_luminal | e_luminal = e_p,i) * P(e_p,i)
p1 <- list(
    crispr = (
        crispr[Description == "FOXA1" & Score < crispr[Line == "T47D", Score], .N]
        / crispr[Description == "FOXA1", .N]
    ),
    rnai = 1/2 * sum(
        rnai[Description == "FOXA1" & Score < rnai[Line == "EFM19", Score], .N],
        rnai[Description == "FOXA1" & Score < rnai[Line == "MCF7", Score], .N]
    ) / rnai[Description == "FOXA1", .N]
)
print(p1)

# Test 2: For m randomly selected lines, what is the probability that
# median({e_i}) < median({e_luminal})

medians <- list(
    crispr = crispr[Description == "FOXA1" & Line == "T47D", median(Score)],
    rnai = rnai[Description == "FOXA1" & Line %in% c("EFM19", "MCF7"), median(Score)]
)

# number of combinations for each method
n_perm <- list(
    # should be equivalent to Test 1
    crispr = choose(crispr[Description == "FOXA1", .N], 1),
    # should be different to Test 1
    rnai = choose(rnai[Description == "FOXA1", .N], 2)
)
cat(paste(n_perm$crispr, "permutations for CRISPR-based essentiality\n"))
cat(paste(n_perm$rnai, "permutations for RNA-based essentiality\n"))

# CRISPR combinations
cat("Performing CRISPR combinations\n")
crispr_combs <- data.table(
    t(combn(crispr[, unique(Line)], 2))
)
colnames(crispr_combs) <- c("Line1", "Line2")
crispr_combs[, Median := mapply(
    function(line1, line2) {
        median(c(
            crispr[Description == "FOXA1" & Line == line1, Score],
            crispr[Description == "FOXA1" & Line == line2, Score]
        ))
    },
    Line1,
    Line2
)]

# RNAi combinations
cat("Performing RNAi combinations\n")
rnai_foxa1 <- rnai[Description == "FOXA1", .(Line, Score)]
rnai_combs <- data.table(
    t(combn(rnai[, unique(Line)], 3))
)
# combine cell line combinations with Scores
rnai_combs <- merge(
	x = rnai_combs,
	y = rnai_foxa1,
	by.x = "V1",
	by.y = "Line"
)
rnai_combs <- merge(
	x = rnai_combs,
	y = rnai_foxa1,
	by.x = "V2",
	by.y = "Line"
)
rnai_combs <- merge(
	x = rnai_combs,
	y = rnai_foxa1,
	by.x = "V3",
	by.y = "Line"
)
colnames(rnai_combs) <- c("Line3", "Line2", "Line1", "Score1", "Score2", "Score3")
rnai_combs_scores <- rnai_combs[, median(c(Score1, Score2, Score3)), by = 1:nrow(rnai_combs)]

p2 <- list(
    crispr = crispr_combs[Median < medians$crispr, .N] / n_perm$crispr,
    rnai = rnai_combs_scores[V1 < medians$rnai, .N] / n_perm$rnai
)

p_vals <- data.table(
    Test = rep(c(1, 2), each = 2),
    Method = rep(c("CRISPR", "RNAi"), 2),
    P = c(unlist(p1), unlist(p2))
)

fwrite(
    p_vals,
    "essentiality-combinations_luminal.tsv",
    col.names = TRUE,
    sep = "\t"
)
