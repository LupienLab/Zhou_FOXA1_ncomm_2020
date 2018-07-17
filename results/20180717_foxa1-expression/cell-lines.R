# ==============================================================================
# Environment
# ==============================================================================
suppressMessages(library("data.table"))
suppressMessages(library("ggplot2"))

# ==============================================================================
# Data
# ==============================================================================
# read CCLE data
ccle <- fread(
    "../../data/external/CCLE_Expression_2012-09-29.res",
    header = TRUE,
    sep = "\t",
    blank.lines.skip = TRUE,  # second line of the res file is blank
    fill = TRUE  # third line of res file only has one column
)

# as per the RES format, the 3rd line contains the number of rows in the file
# remove this line from the table
ccle <- ccle[-1]

# as per the RES format, every second column is the Absent/Present/Marginal call
# for the previous column. Ignore these columns and remove them
# (these columns are all preface with 'V' since fread assigns them a dummy name)
cell_line_idx <- grepl("^V\\d+", colnames(ccle))
ccle <- ccle[, which(cell_line_idx) := NULL]

# select only expression columns relating to prostate cell lines and
# melt data into 2 columns, one for the cell line, one for the expression value
prostate_idx <- grep("PROSTATE", colnames(ccle))
prostate_idx <- c(1, 2, prostate_idx)
pca_exprs <- melt(
    ccle[, ..prostate_idx],  # reduce to only columns of interest
    id.vars = c(1, 2)  # melt based on first two columns
)
colnames(pca_exprs) <- c("Description", "Accession", "Cell", "Expression")

# ==============================================================================
# Analysis
# ==============================================================================
# ==============================================================================
# Plots
# ==============================================================================
gg <- (
    ggplot(data = pca_exprs, aes(x = Cell, y = Expression, fill = Cell))
    + geom_boxplot()
    + labs(x = "Cell Line", y = "Gene Expression (RMA-normalized)")
)
ggsave(
    "cell-lines.png",
    height = 25,
    width = 40,
    units = "cm",
    bg = "transparent"
)
