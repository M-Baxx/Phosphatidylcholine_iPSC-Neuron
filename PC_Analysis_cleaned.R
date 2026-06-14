# ============================================================================
#  Phosphatidylcholine (PC) lipidomics analysis
#  iPSC  vs  neural-induced (NI)  vs  differentiated (DA)
#  across three iPSC lines: HC1, HC2, SAD
#
#  Pipeline:
#    load -> filter -> log2 -> impute -> per-group mean/SD -> normalise
#         -> pairwise t-tests (volcano tables) -> correlation/heatmap -> PCA
#
#  Input data: Supplementary Data accompanying the manuscript.
#    Expected: a CSV whose first column holds PC species names, plus one column
#    per sample. Sample columns are matched BY NAME (see `samples` below), so
#    column order / extra columns in the file do not matter.
#
#  All outputs (tables + figures) are written to `output_dir`.
# ============================================================================

# ---- Configuration ---------------------------------------------------------
data_dir   <- "data"                                   # folder with input CSV
output_dir <- "output"                                 # generated tables/figures
input_csv  <- file.path(data_dir, "PC_iPSC_NI_DA.csv") # <- supplementary file
write_intermediate <- TRUE     # FALSE = write only key outputs, skip per-group CSVs

# NOTE (reproducibility): imputation draws random replacements for missing
# values. The original analysis was NOT seeded ("if you redo this you will get
# different newly-assigned values"), so the exact published imputed values
# cannot be regenerated. A seed is set here so this script is reproducible
# going forward; treat the original imputed CSV as the reference for the
# figures already in the manuscript.
set.seed(2020)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(output_dir, f)            # helper for output paths

# ---- Packages --------------------------------------------------------------
# install.packages(c("ggplot2", "gplots"))            # uncomment on first run
library(ggplot2)   # volcano + PCA plots
library(gplots)    # heatmap.2
# (Original also loaded VennDiagram, reshape2 and plyr; none are used here.)

# ---- Sample sheet (single source of truth) ---------------------------------
# 38 samples in canonical order. Everything downstream selects columns BY NAME
# from this sheet rather than by position, removing the index fragility of the
# original (e.g. subset(..., c(6:8))).
sample_counts <- list(
  HC1 = c(iPSC = 5, NI = 3, DA = 4),
  HC2 = c(iPSC = 4, NI = 5, DA = 4),
  SAD = c(iPSC = 4, NI = 5, DA = 4)
)
states <- c("iPSC", "NI", "DA")

samples <- do.call(rbind, lapply(names(sample_counts), function(line) {
  do.call(rbind, lapply(states, function(state) {
    n <- sample_counts[[line]][[state]]
    data.frame(line = line, state = state, rep = seq_len(n),
               name = paste(line, state, seq_len(n)), stringsAsFactors = FALSE)
  }))
}))
# samples$name -> "HC1 iPSC 1" ... "SAD DA 4"  (38 rows)

# ---- Load input ------------------------------------------------------------
raw <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
species <- raw[[1]]                                    # first column = PC species
# tolerate "HC1.iPSC.1" / "HC1_iPSC_1" style headers by normalising to spaces
colnames(raw) <- gsub("[._]+", " ", trimws(colnames(raw)))

missing <- setdiff(samples$name, colnames(raw))
if (length(missing))
  stop("Input is missing expected sample columns:\n  ", paste(missing, collapse = "\n  "))

intens <- as.matrix(raw[, samples$name, drop = FALSE])  # 38 samples, canonical order

# small helpers to pull a group's replicate columns by name
cols_of  <- function(line, state) samples$name[samples$line == line & samples$state == state]
reps_of  <- function(line, state) intens_df[, cols_of(line, state), drop = FALSE]
state_reps <- function(state)     intens_df[, samples$name[samples$state == state], drop = FALSE]

# ============================================================================
#  1. Filter
#     Within each (line, state) group, blank a species (set its values to 0)
#     if it has >= 2 zero values among that group's replicates.
# ============================================================================
filtered <- intens
for (line in names(sample_counts)) for (state in states) {
  cc   <- cols_of(line, state)
  g    <- intens[, cc, drop = FALSE]
  keep <- rowSums(g == 0) < 2                          # FALSE -> row blanked to 0
  filtered[, cc] <- sweep(g, 1, keep, `*`)
}

# ============================================================================
#  2. log2 transform  (zeros -> -Inf -> NA)
# ============================================================================
logmat <- log2(filtered)
logmat[!is.finite(logmat)] <- NA

# ============================================================================
#  3. Impute missing values
#     For each sample column, replace NAs with draws from
#     N(min(col), sd(lower quartile of col)).
#
#  NOTE: the original imputed object columns 1:39, where column 1 was the
#  (non-numeric) species column and used quantile()[1]; the 38 real data
#  columns (2:39) all used quantile()[2] (the 25th percentile). This loop
#  reproduces the data-affecting behaviour by imputing the 38 numeric sample
#  columns uniformly with quantile()[2].
# ============================================================================
imputed <- as.data.frame(logmat)
for (j in seq_along(imputed)) {
  na_idx <- is.na(imputed[[j]])
  if (any(na_idx)) {
    lowq <- imputed[[j]] <= quantile(imputed[[j]], na.rm = TRUE)[2]   # <= 25th pct
    imputed[na_idx, j] <- rnorm(sum(na_idx),
                                mean = min(imputed[[j]], na.rm = TRUE),
                                sd   = sd(imputed[lowq, j], na.rm = TRUE))
  }
}
intens_df <- imputed                                   # used by reps_of()/state_reps()
rownames(intens_df) <- species
if (write_intermediate) write.csv(intens_df, out("Lipids_Imputed.csv"))

# ============================================================================
#  4. Per-group means and SDs (optional intermediate tables)
# ============================================================================
if (write_intermediate) {
  for (line in names(sample_counts)) for (state in states) {
    g <- reps_of(line, state)
    write.csv(cbind(g, mean = rowMeans(g)),   out(sprintf("%s_%s_mean.csv", line, state)))
    write.csv(cbind(g, sd   = apply(g, 1, sd)), out(sprintf("%s_%s_sd.csv", line, state)))
  }
}

# ============================================================================
#  5. Normalise
#     Centre each replicate on the mean of the OTHER two states (same line),
#     i.e. scaled = replicate - mean(other_state_A_mean, other_state_B_mean).
#     Output layout per group: [ replicates | mean | Scaled.1 ... Scaled.n ].
# ============================================================================
make_normalised <- function(reps, other_a_mean, other_b_mean) {
  centre <- (other_a_mean + other_b_mean) / 2
  scaled <- as.data.frame(lapply(reps, function(col) col - centre))
  names(scaled) <- paste0("Scaled.", seq_along(reps))
  cbind(reps, mean = rowMeans(reps), scaled)
}

normalised <- list()
# per cell line
for (line in names(sample_counts)) {
  m <- setNames(lapply(states, function(s) rowMeans(reps_of(line, s))), states)
  for (s in states) {
    o <- setdiff(states, s)
    normalised[[paste(line, s)]] <- make_normalised(reps_of(line, s), m[[o[1]]], m[[o[2]]])
  }
}
# pooled across lines ("All")
allm <- setNames(lapply(states, function(s) rowMeans(state_reps(s))), states)
for (s in states) {
  o <- setdiff(states, s)
  normalised[[paste("All", s)]] <- make_normalised(state_reps(s), allm[[o[1]]], allm[[o[2]]])
}

# ============================================================================
#  6. Pairwise t-tests / volcano tables
#     Welch t-test per species, BH-adjusted, with fold change and thresholds.
#     Thresholds: |log2 FC| >= 1.5 and -log10(adj p) > 1.3.
# ============================================================================
scaled_of <- function(key) {                           # Scaled.* columns of a group
  df <- normalised[[key]]
  df[, grep("^Scaled\\.", names(df)), drop = FALSE]
}

volcano_table <- function(a, b, lfc = 1.5, neglog10p = 1.3) {
  na_ <- ncol(a); nb_ <- ncol(b); D <- cbind(a, b)
  D$P.value <- apply(D, 1, function(r)
    t.test(x = r[1:na_], y = r[(na_ + 1):(na_ + nb_)],
           alternative = "two.sided", var.equal = FALSE)$p.value)
  D$adj.pvalues          <- p.adjust(D$P.value, method = "BH", n = length(D$P.value))
  D$normalised.mean.A    <- apply(D[1:na_], 1, mean)
  D$Abundance.A          <- 2 ^ D$normalised.mean.A
  D$normalised.mean.B    <- apply(D[(na_ + 1):(na_ + nb_)], 1, mean)
  D$Abundance.B          <- 2 ^ D$normalised.mean.B
  D$Fold.change          <- D$Abundance.A / D$Abundance.B
  D$Log.Fold.change      <- log2(D$Fold.change)
  D$neg.log10.adj.pvalue <- -log10(D$adj.pvalues)
  D$Threshold     <- as.factor(D$Log.Fold.change >=  lfc & D$neg.log10.adj.pvalue > neglog10p)
  D$Threshold.neg <- as.factor(D$Log.Fold.change <= -lfc & D$neg.log10.adj.pvalue > neglog10p)
  D
}

# group A (first) vs group B (second) for each comparison
comparisons <- list(
  c("HC1 iPSC", "HC1 NI"), c("HC1 iPSC", "HC1 DA"), c("HC1 NI", "HC1 DA"),
  c("HC2 iPSC", "HC2 NI"), c("HC2 iPSC", "HC2 DA"), c("HC2 NI", "HC2 DA"),
  c("SAD iPSC", "SAD NI"), c("SAD iPSC", "SAD DA"), c("SAD NI", "SAD DA"),
  c("All iPSC", "All NI"), c("All iPSC", "All DA"), c("All NI", "All DA")
)

volcano <- list()
for (cmp in comparisons) {
  key <- paste0(sub(" ", ".", cmp[1]), ".vs.", sub("^\\S+ ", "", cmp[2]))
  vt  <- volcano_table(scaled_of(cmp[1]), scaled_of(cmp[2]))
  rownames(vt) <- species
  volcano[[key]] <- vt
  write.csv(vt, out(paste0("Ttest_", key, ".csv")))
}

# ---------------------------------------------------------------------------
#  DISCREPANCIES vs the original RAW script (please reconcile with the paper)
#  The generic volcano_table() above is internally consistent. The original
#  hand-written blocks contained a few copy-paste slips that this version does
#  NOT reproduce; if the manuscript figures used the original behaviour, revert
#  the affected comparison(s) by hand:
#    - iPSC vs DA (pooled "All"): original set Abundance.iPSC = 2^mean.DA
#      (referenced the wrong mean). Here it correctly uses mean.iPSC.
#    - NI vs DA (pooled "All"):   original set Abundance.NI = 2^mean.DA, and
#      assigned $Threshold twice (the negative test overwrote the positive one).
#      Here both abundances and both thresholds are computed correctly.
#    - HC1 NI vs DA: original Threshold.neg used "<= 1.5" instead of "<= -1.5".
#      Here it uses -1.5, consistent with every other comparison.
# ---------------------------------------------------------------------------

# ============================================================================
#  7. Assembled normalised matrix (38 samples) for correlation / heatmap / PCA
#     Built from the per-line scaled values, in canonical sample order.
# ============================================================================
lipids_normalised <- do.call(cbind, lapply(seq_len(nrow(samples)), function(i) {
  scaled_of(paste(samples$line[i], samples$state[i]))[, samples$rep[i], drop = FALSE]
}))
colnames(lipids_normalised) <- samples$name
rownames(lipids_normalised) <- species
write.csv(lipids_normalised, out("Lipids_Normalised.csv"))

mat_norm <- as.matrix(lipids_normalised)

# ---- Pearson correlation between samples ----
all_pearson <- cor(mat_norm, method = "pearson", use = "complete.obs")
write.csv(all_pearson, out("All_Pearson.csv"))

# ---- Heatmap with dendrograms ----
# NOTE: the original passed `dendogram=` (misspelled, so silently ignored) and
# cexRow = 0.01 (row labels invisible). Corrected here.
pdf(out("Heatmap_Lipids_Normalised.pdf"), width = 8, height = 10)
heatmap.2(data.matrix(lipids_normalised), main = "", notecol = "black",
          density.info = "none", trace = "none", col = "bluered",
          dendrogram = "both", srtCol = 45, cexCol = 0.75, cexRow = 0.4, keysize = 1)
dev.off()

# ============================================================================
#  8. PCA
#     Lipids = rows, samples = columns; samples are positioned using the
#     rotation (loadings), as in the original. Axes show % variance; states are
#     relabelled NI -> NR and DA -> N; 95% CI ellipses are drawn per state.
# ============================================================================
# point colours by line x state (matches the original palette hues)
point_palette <- c(
  "HC1 iPSC" = "#b69dff", "HC2 iPSC" = "#8f8dff", "SAD iPSC" = "#5762ff",
  "HC1 NR"   = "#9dd5ff", "HC2 NR"   = "#8dfffd", "SAD NR"   = "#0a6699",
  "HC1 N"    = "#b5ffb5", "HC2 N"    = "#5fff57", "SAD N"    = "#1f9e1f"
)
state_fills <- c("iPSC" = "#8f8dff", "NR" = "#57a8ff", "N" = "#5fff57")

pca_plot <- function(mat, file_prefix, drop_lines = NULL) {
  if (!is.null(drop_lines)) {
    pat <- paste0("^(", paste(drop_lines, collapse = "|"), ")")
    mat <- mat[, !grepl(pat, colnames(mat)), drop = FALSE]
  }
  mat <- mat[apply(mat, 1, sd) > 0, ]                  # drop zero-variance lipids
  pc  <- prcomp(mat, center = TRUE, scale. = TRUE)
  pv  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)

  L   <- data.frame(pc$rotation)
  lab <- rownames(L)
  lab <- gsub(" NI ", " NR ", lab)
  lab <- gsub(" DA ", " N ",  lab)

  line  <- sub(" .*$", "", lab)
  L$State <- factor(sub("^[^ ]+ ([^ ]+).*$", "\\1", lab), levels = c("iPSC", "NR", "N"))
  point_cols <- point_palette[paste(line, L$State)]

  disp <- lab
  disp <- sub("^HC1 ", "Line 1 ", disp)
  disp <- sub("^HC2 ", "Line 2 ", disp)
  disp <- sub("^SAD ", "Line 3 ", disp)

  g <- ggplot(L, aes(PC1, PC2)) +
    stat_ellipse(aes(fill = State), geom = "polygon",
                 alpha = 0.15, colour = NA, type = "t", level = 0.95) +
    geom_point(size = 5, colour = point_cols) +
    geom_text(aes(label = disp), hjust = 0, vjust = 2, size = 3) +
    scale_fill_manual(values = state_fills) +
    labs(title = "",
         x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)"),
         fill = "State") +
    theme(panel.background = element_rect(fill = "white"),
          panel.border     = element_blank(),
          axis.line        = element_line(colour = "black", linewidth = 0.8),
          plot.title       = element_text(face = "bold", size = rel(2)))

  ggsave(out(paste0(file_prefix, ".pdf")), g, width = 7, height = 6)
  ggsave(out(paste0(file_prefix, ".png")), g, width = 7, height = 6, dpi = 600)
  g
}

pca_plot(mat_norm, "PCA_Lipids_Normalised")                       # all samples
pca_plot(mat_norm, "PCA_Lipids_Lines1-2_noSAD", drop_lines = "SAD")  # HC1 + HC2 only

# ============================================================================
#  End of pipeline. All tables and figures are in `output_dir`.
# ============================================================================
