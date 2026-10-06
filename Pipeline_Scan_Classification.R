# ======================================================================
# PIPELINE SCAN CLASSIFICATION
# Binary classification of iron and copper in wine relative to an
# enological threshold (10 and 1 mg/L, respectively) from FT-MIR spectra.
# See the README for a complete description of the workflow,
# the input data and the generated outputs.
#
# Models evaluated: LDA, PLS-DA, SVM, RF, CART, kNN, XGBoost, Naive Bayes
# Metrics: Accuracy, Kappa, Sensitivity, Specificity, F1, AUC
# ======================================================================

# ======================================
# PREPROCESSING CONFIGURATION
# ======================================
# SCATTER_FIRST: order of operations, always applied to the FULL spectrum.
#   TRUE  -> scatter correction (SNV/MSC) -> Savitzky-Golay (smoothing or derivative)
#            SNV and MSC were designed for non-derivatized absorbance spectra.
#   FALSE -> Savitzky-Golay -> SNV/MSC (alternative order).
# SMOOTH_ALWAYS: if TRUE, combinations WITHOUT a derivative are also smoothed with
#   Savitzky-Golay (derivative order 0, polynomial order 3, 11-point window). If FALSE,
#   they are left as unsmoothed spectra.
SCATTER_FIRST <- TRUE
SMOOTH_ALWAYS <- TRUE

# Name of each combination (tables, plots and files). SG0 = smoothing only;
# SG1 / SG2 = 1st / 2nd derivative. Techniques appear in the name in the order in
# which they were applied (e.g. "SNV + SG1" = SNV followed by 1st derivative).
make_pp_name <- function(sc, dv) {
  sg      <- if (dv == 0 && !SMOOTH_ALWAYS) NULL else paste0("SG", dv)
  sc_part <- if (sc == "none") NULL else sc
  parts   <- if (SCATTER_FIRST) c(sc_part, sg) else c(sg, sc_part)
  if (length(parts) == 0) "raw" else paste(parts, collapse = " + ")
}

# ======================================
# 0. Parallelization
# ======================================
library(doParallel)
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

cat("\n", rep("=", 70), "\n", sep="")
cat("FT-MIR BINARY CLASSIFICATION PIPELINE (Scan Classification)\n")
cat(rep("=", 70), "\n", sep="")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("CPU cores:", num_cores, "\n")
cat(rep("=", 70), "\n\n", sep="")

# ======================================
# 1. Libraries
# ======================================
library(dplyr)
library(tidyr)
library(caret)
library(pracma)
library(prospectr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(scales)
library(pROC)
library(Boruta)
library(MASS)
library(pls)
library(e1071)
library(randomForest)
library(rpart)
library(xgboost)
library(naivebayes)
library(Matrix)
library(dplyr)

library(dplyr)
library(MASS)

select <- dplyr::select
filter <- dplyr::filter
lag <- dplyr::lag

# Resolve name conflict for the margin function
margin <- ggplot2::margin

# ======================================
# GLOBAL SETTINGS
# ======================================
POS_CLASS <- NULL  # assigned automatically when the data are loaded

# ======================================
# 2. Directories
# ======================================
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

for (d in c("Results", "Results/Excel", "Results/Spectra",
            "Results/Heatmaps", "Results/Other",
            "Results/Outliers", "Results/ConfusionMatrices",
            "Results/ROC", "Results/ScatterPlots")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
cat(">>> Output folders created.\n\n")

# ======================================
# 3. Data loading
# ======================================
cat(">>> Loading data...\n")

data_path  <- "data/FINAL_DATA_SET.xlsx"      # path to the data file (see the README)
sheet_name <- "IRON"   # analyte to run: "IRON" or "COPPER"

# Locates the data file regardless of the working directory or of whether the
# name uses spaces or underscores ("FINAL DATA SET.xlsx" / "FINAL_DATA_SET.xlsx").
locate_data_file <- function(preferred) {
  if (file.exists(preferred)) return(preferred)
  found <- list.files(".", pattern = "^FINAL[ _]DATA[ _]SET\\.xlsx$", recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  if (length(found) > 0) return(found[1])
  stop("ERROR: data file not found. Current directory: ", getwd(),
       "\n  Copy 'FINAL_DATA_SET.xlsx' to that folder (or to 'data/') or adjust data_path.")
}
data_path <- locate_data_file(data_path)

# Accepts sheet names in English or Spanish (and the 'COOPER' variant found in the dataset).
resolve_sheet <- function(path, wanted) {
  syn <- list(POTASSIUM = c("POTASSIUM", "POTASIO"), MAGNESIUM = c("MAGNESIUM", "MAGNESIO"),
              CALCIUM   = c("CALCIUM", "CALCIO"),    IRON      = c("IRON", "HIERRO"),
              COPPER    = c("COPPER", "COOPER", "COBRE"))
  have <- getSheetNames(path)
  key  <- names(syn)[vapply(syn, function(v) toupper(wanted) %in% v, logical(1))]
  cand <- if (length(key) == 1) syn[[key]] else toupper(wanted)
  hit  <- have[toupper(have) %in% cand]
  if (length(hit) == 0)
    stop("ERROR: sheet not found for '", wanted, "'. Available sheets: ", paste(have, collapse = ", "))
  hit[1]
}
sheet_in_file <- resolve_sheet(data_path, sheet_name)
cat("   Data file:", data_path, "| sheet:", sheet_in_file, "\n")

nir        <- read.xlsx(data_path, sheet = sheet_in_file)
nir        <- data.frame(lapply(nir, function(x) if (is.character(x)) as.factor(x) else x))
sample_ids <- as.character(nir[, 1])

y_raw      <- as.factor(nir[, 2])
levels(y_raw) <- make.names(levels(y_raw))
y_raw      <- factor(y_raw)

# The POSITIVE class is the one of interest: samples ABOVE the threshold ("Higher...").
# The whole pipeline (PLS-DA threshold, XGBoost 0/1 coding, roc(levels = ...))
# assumes the positive class is the 2nd factor level, so the levels are reordered
# so that "Higher" comes second (alphabetical order would make "Lower" the
# positive class and invert the meaning of sensitivity, specificity and F1).
lv_all <- levels(y_raw)
pos_lv <- lv_all[grepl("^Higher", lv_all)]
if (length(pos_lv) == 1) {
  y_raw <- factor(y_raw, levels = c(setdiff(lv_all, pos_lv), pos_lv))
} else {
  warning("Class 'Higher...' not identified: using alphabetical order (2nd level = positive).")
}
class_names <- levels(y_raw)

if (is.null(POS_CLASS)) POS_CLASS <- class_names[2]

X_raw       <- nir[, -c(1, 2)]                     # FULL, continuous spectrum

# --------------------------------------------------------------------
# Helper: parses wavenumbers from column names such as "X960.45"
# or "X960,45" (decimal comma, regional Excel format). Without this
# replacement, as.numeric("960,45") silently returns NA and
# breaks patz_idx, the post-preprocessing trimming and the plots with
# Boruta lines. It is used throughout the pipeline instead of a bare
# as.numeric(gsub("X","",...)).
# --------------------------------------------------------------------
parse_wn <- function(x) as.numeric(gsub(",", ".", gsub("X", "", x), fixed = TRUE))

# --------------------------------------------------------------------
# Helpers: format the optimized hyperparameters of each winning model
# (caret$bestTune, or the grid selected by manual CV in XGB)
# as a readable string, to report them together with the performance
# metrics of each combination.
# --------------------------------------------------------------------
format_bestTune <- function(bt) {
  if (is.null(bt) || nrow(bt) == 0) return(NA_character_)
  vals <- sapply(bt, function(x) if (is.numeric(x)) format(round(x, 5), trim = TRUE) else as.character(x))
  paste(paste0(names(bt), "=", vals), collapse = "; ")
}
format_xgb_params <- function(best_params) {
  if (is.null(best_params)) return(NA_character_)
  p <- best_params$params
  paste0("nrounds=", best_params$nrounds,
         "; max_depth=", p$max_depth, "; eta=", p$eta, "; gamma=", p$gamma,
         "; subsample=", p$subsample, "; colsample_bytree=", p$colsample_bytree,
         "; min_child_weight=", p$min_child_weight)
}

wavelengths <- parse_wn(colnames(X_raw))
if (anyNA(wavelengths)) {
  cat("   WARNING:", sum(is.na(wavelengths)),
      "column name(s) could not be parsed as wavenumbers (check decimal separator).\n")
}

cat("   Samples (raw):", nrow(X_raw), "\n")
cat("   Spectral variables (full spectrum):", ncol(X_raw), "\n")
cat("   Classes:", paste(class_names, collapse = " vs "), "\n")
cat("   Positive class (AUC/ROC):", POS_CLASS, "\n\n")

# --------------------------------------------------------------------
# Spectral windows of Patz et al. (2004), retained for modelling.
# The spectrum is trimmed AFTER smoothing/derivatization to avoid artifacts
# at the junctions of non-contiguous windows.
# --------------------------------------------------------------------
PATZ_WINDOWS <- list(c(965, 1582), c(1698, 2006), c(2701, 2971))
in_patz_windows <- function(w) {
  Reduce(`|`, lapply(PATZ_WINDOWS, function(rng) w >= rng[1] & w <= rng[2]))
}
patz_idx <- in_patz_windows(wavelengths)
cat("   Spectral variables within Patz windows:", sum(patz_idx), "\n\n")

# ======================================
# 4. Outlier detection and removal
#    AND criterion: T2 AND Q  (on the Patz-trimmed variables)
# ======================================
cat(">>> Outlier detection (T2 AND Q)...\n")

X_raw_patz <- X_raw[, patz_idx, drop = FALSE]
X_scaled <- scale(X_raw_patz)
pca_exp  <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
var_cum  <- cumsum(pca_exp$sdev^2) / sum(pca_exp$sdev^2)
n_pcs    <- max(3, which(var_cum >= 0.95)[1])
cat("   PCs used:", n_pcs, sprintf("(%.1f%% variance)\n", var_cum[n_pcs] * 100))

scores      <- pca_exp$x[, 1:n_pcs, drop = FALSE]
var_exp2    <- round(pca_exp$sdev^2 / sum(pca_exp$sdev^2) * 100, 1)

lambda_inv  <- diag(1 / pca_exp$sdev[1:n_pcs]^2)
T2          <- rowSums((scores %*% lambda_inv) * scores)
T2_lim      <- qchisq(0.99, df = n_pcs)
flag_T2     <- T2 > T2_lim

X_rec       <- scores %*% t(pca_exp$rotation[, 1:n_pcs, drop = FALSE])
Q           <- rowSums((X_scaled - X_rec)^2)
Q_lim       <- mean(Q) + 3 * sd(Q)
flag_Q      <- Q > Q_lim

flag_outlier <- flag_T2 & flag_Q   # AND: only "severe" outliers that violate BOTH criteria at once
n_outliers   <- sum(flag_outlier)

cat("   T2 flagged:", sum(flag_T2), "| Q flagged:", sum(flag_Q),
    "| Removed (T2 AND Q):", n_outliers, "of", nrow(X_raw), "\n")

if (n_outliers > 0) {
  oidx    <- which(flag_outlier)
  out_rep <- data.frame(
    Sample_ID  = sample_ids[oidx], Row_Index = oidx,
    Class      = as.character(y_raw[oidx]),
    T2_score   = round(T2[oidx], 3),  T2_limit  = round(T2_lim, 3), T2_flagged = flag_T2[oidx],
    Q_residual = round(Q[oidx], 4),   Q_limit   = round(Q_lim, 4),  Q_flagged  = flag_Q[oidx],
    Criterion  = "T2 > chi2(0.99) AND Q > mean+3SD", stringsAsFactors = FALSE)
  write.xlsx(out_rep, "Results/Outliers/Outliers_Removed.xlsx", overwrite = TRUE)
  cat("   Removed:", paste(out_rep$Sample_ID, collapse = ", "), "\n")
} else {
  write.xlsx(data.frame(Message = "No outliers detected.", T2_limit = round(T2_lim, 3),
                        Q_limit = round(Q_lim, 4), Criterion = "T2 AND Q"),
             "Results/Outliers/Outliers_Removed.xlsx", overwrite = TRUE)
  cat("   No outliers detected.\n")
}

# Color palette by class
cls_colors <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                          "#FF7F00", "#A65628")[seq_len(length(class_names))],
                        class_names)
# Fixed colors by meaning (independent of level order):
# red = above the threshold (positive class), blue = below.
cls_colors[grepl("^Higher", class_names)] <- "#E41A1C"
cls_colors[grepl("^Lower",  class_names)] <- "#377EB8"

# Raw spectra plot
cat("   Plotting raw spectra...\n")
raw_long <- data.frame(
  Wavenumber = rep(wavelengths, times = nrow(X_raw)),
  Absorbance = as.vector(t(as.matrix(X_raw))),
  Sample_ID  = rep(sample_ids, each = length(wavelengths)),
  Class      = rep(as.character(y_raw), each = length(wavelengths)),
  Outlier    = rep(ifelse(flag_outlier, "Outlier", "Normal"), each = length(wavelengths)),
  stringsAsFactors = FALSE)

raw_normal  <- subset(raw_long, Outlier == "Normal")
raw_outlier <- subset(raw_long, Outlier == "Outlier")

p_raw <- ggplot() +
  geom_line(data = raw_normal,
            aes(x = Wavenumber, y = Absorbance, group = Sample_ID, color = Class),
            alpha = 0.35, linewidth = 0.4) +
  geom_line(data = raw_outlier,
            aes(x = Wavenumber, y = Absorbance, group = Sample_ID),
            color = "black", alpha = 0.85, linewidth = 0.7, linetype = "dashed") +
  scale_color_manual(values = cls_colors, name = "Class") +
  labs(title    = "Raw FTMIR Spectra",
       subtitle = paste0("All samples before preprocessing  |  ",
                         paste(paste0(class_names, " (n=",
                                      table(y_raw[!flag_outlier]), ")"),
                               collapse = "  |  "),
                         "  |  Outliers removed (dashed, n=", n_outliers, ")"),
       x = expression("Wavenumber (cm"^{-1}*")"), y = "Absorbance") +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
        panel.grid.minor = element_blank(),
        panel.border  = element_rect(color = "gray70", fill = NA))

if (n_outliers > 0) {
  olbl <- raw_outlier %>% group_by(Sample_ID) %>% slice_max(Absorbance, n = 1) %>% ungroup()
  p_raw <- p_raw +
    geom_label_repel(data = olbl, aes(x = Wavenumber, y = Absorbance, label = Sample_ID),
                     color = "black", fill = "white", fontface = "bold",
                     size = 3, box.padding = 0.4, max.overlaps = 20)
}
ggsave("Results/Outliers/Raw_Spectra_Outliers.png", p_raw, width = 12, height = 6, dpi = 300, bg = "white")
cat("   OK Raw spectra saved\n")

# PCA plots
pca_df <- data.frame(
  PC1 = pca_exp$x[, 1], PC2 = pca_exp$x[, 2], PC3 = pca_exp$x[, 3],
  Sample_ID = sample_ids, Class = as.character(y_raw),
  Outlier   = ifelse(flag_outlier, "Outlier", "Normal"), stringsAsFactors = FALSE)

make_pca_cls <- function(df, xc, yc, xlab, ylab, title, fname) {
  dn <- subset(df, Outlier == "Normal")
  do <- subset(df, Outlier == "Outlier")
  p <- ggplot(df, aes(x = .data[[xc]], y = .data[[yc]])) +
    geom_point(data = dn, aes(color = Class), size = 3.5, alpha = 0.85) +
    scale_color_manual(values = cls_colors, name = "Class") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          legend.position = "right", panel.border = element_rect(color = "gray70", fill = NA))
  if (nrow(do) > 0)
    p <- p +
      geom_point(data = do, color = "black", fill = "red", shape = 23, size = 4.5, stroke = 1.2) +
      annotate("text", x = min(df[[xc]], na.rm = TRUE), y = max(df[[yc]], na.rm = TRUE),
               label = paste0("Outliers removed: ", nrow(do)),
               hjust = 0, vjust = 1, color = "red", fontface = "bold", size = 3.8)
  ggsave(fname, p, width = 11, height = 8, dpi = 300, bg = "white")
  cat("   OK PCA saved:", basename(fname), "\n")
  invisible(p)
}

make_pca_cls(pca_df, "PC1", "PC2",
             paste0("PC1 (", var_exp2[1], "%)"), paste0("PC2 (", var_exp2[2], "%)"),
             "Exploratory PCA — PC1 vs PC2", "Results/Outliers/PCA_PC1_vs_PC2.png")
make_pca_cls(pca_df, "PC2", "PC3",
             paste0("PC2 (", var_exp2[2], "%)"), paste0("PC3 (", var_exp2[3], "%)"),
             "Exploratory PCA — PC2 vs PC3", "Results/Outliers/PCA_PC2_vs_PC3.png")

# Influence plot
infl_df <- data.frame(T2 = T2, Q = Q, Sample_ID = sample_ids,
                       Outlier = ifelse(flag_outlier, "Outlier", "Normal"), stringsAsFactors = FALSE)
p_infl <- ggplot(infl_df, aes(x = T2, y = Q, color = Outlier)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(data = subset(infl_df, Outlier == "Normal"), aes(label = Sample_ID),
                  size = 3.4, color = "gray40", box.padding = 0.2, max.overlaps = 20) +
  geom_vline(xintercept = T2_lim, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_hline(yintercept = Q_lim,  linetype = "dashed", color = "steelblue", linewidth = 0.8) +
  scale_color_manual(values = c("Normal" = "steelblue", "Outlier" = "red")) +
  labs(title    = "Influence Plot: Hotelling T2 vs Q Residuals",
       subtitle = "Outliers: T2 > chi2(0.99) AND Q > mean+3SD",
       x = "Hotelling T2", y = "Q Residuals", color = "") +
  theme_minimal(base_size = 22) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 24),
        plot.subtitle = element_text(hjust = 0.5, size = 17, color = "gray40"),
        axis.title = element_text(size = 20), axis.text = element_text(size = 17),
        legend.text = element_text(size = 17), legend.title = element_text(size = 18),
        legend.position = "top", panel.border = element_rect(color = "gray70", fill = NA))
if (n_outliers > 0)
  p_infl <- p_infl +
    geom_label_repel(data = subset(infl_df, Outlier == "Outlier"), aes(label = Sample_ID),
                     color = "red", fill = "white", fontface = "bold",
                     size = 4.5, box.padding = 0.4, max.overlaps = 30)
ggsave("Results/Outliers/Influence_Plot_T2_vs_Q.png", p_infl, width = 11, height = 8, dpi = 300, bg = "white")
cat("   OK Influence plot saved\n")

# ======================================
# 4bis. Circular dendrogram of samples
#   Hierarchical clustering of the spectra (outliers excluded), with the
#   sample ID colored by OIV compliance class (the same binary class
#   used for classification: above or below the critical limit
#   of the corresponding metal).
# ======================================
cat(">>> Building circular dendrogram (samples colored by OIV class)...\n")
if (!requireNamespace("dendextend", quietly = TRUE)) install.packages("dendextend")
if (!requireNamespace("circlize",   quietly = TRUE)) install.packages("circlize")
suppressPackageStartupMessages({library(dendextend); library(circlize)})

X_clean   <- X_raw_patz[!flag_outlier, , drop = FALSE]
y_clean   <- droplevels(y_raw[!flag_outlier])
ids_clean <- sample_ids[!flag_outlier]
rownames(X_clean) <- make.unique(ids_clean)

d_clust  <- dist(scale(X_clean), method = "euclidean")
hc_clust <- hclust(d_clust, method = "ward.D2")
dend     <- as.dendrogram(hc_clust)

# Leaf labels WITHOUT the actual sample ID: the class the sample belongs to
# plus a sequential number within that class
# (e.g. "Higher_01", "Higher_02", ..., "Lower_01", ...), following the style of the
# reference dendrogram (Country_Number) but without exposing the internal
# sample code. Computed BEFORE renaming the leaves,
# while "labels(dend)" still holds the actual rownames, so that the
# match against y_clean/rownames(X_clean) can be made.
leaf_ids    <- labels(dend)
leaf_class  <- y_clean[match(leaf_ids, rownames(X_clean))]
leaf_colors <- cls_colors[as.character(leaf_class)]

class_short  <- setNames(sub("^([A-Za-z]+).*", "\\1", gsub("\\.", " ", class_names)), class_names)
leaf_short   <- class_short[as.character(leaf_class)]
leaf_counter <- ave(seq_along(leaf_short), leaf_short, FUN = seq_along)

# LABEL color = actual class (known ground truth, the same
# cls_colors palette used throughout the pipeline).
labels_colors(dend) <- leaf_colors
labels(dend)         <- sprintf("%s_%02d", leaf_short, leaf_counter)

# BRANCH color = hierarchical cluster found by Ward.D2 (unsupervised),
# to visually compare whether the dendrogram clusters
# match the actual class (label color).
# k_clusters is adjustable: larger k = more distinguished colors/branches.
k_clusters <- 4
dend <- color_branches(dend, k = k_clusters)
dend <- set(dend, "branches_lwd", 2.2)
dend <- set(dend, "labels_cex", 1.05)

png("Results/Outliers/Dendrogram_Circular.png", width = 3600, height = 3600, res = 300)
par(mar = c(1, 1, 5, 1), xpd = TRUE)   # top margin so the title is not cut off
circlize_dendrogram(dend, labels_track_height = 0.28, dend_track_height = 0.55)
title(main = paste0("Hierarchical Clustering of Samples (Ward.D2, Euclidean)\n",
                     "Branch color = ", k_clusters, " hierarchical clusters  |  ",
                     "Label color = actual class"), cex.main = 1.7, line = 1)
legend("bottomright", legend = gsub("\\.", " ", names(cls_colors)), text.col = cls_colors,
       bty = "n", cex = 1.6, pt.cex = 0)
dev.off()
cat("   OK Circular dendrogram saved: Results/Outliers/Dendrogram_Circular.png\n")

# Remove outliers
keep_idx   <- !flag_outlier
X          <- X_raw[keep_idx, , drop = FALSE]
y          <- droplevels(y_raw[keep_idx])
sample_ids <- sample_ids[keep_idx]
cat(sprintf("\n   After removal: %d samples (removed %d)\n\n", sum(keep_idx), n_outliers))

# ======================================
# 5. Train / Test Split (stratified)
# ======================================
# Stratified random sampling (createDataPartition) with a fixed seed.
# Distance-based algorithms (Kennard-Stone, Duplex) are avoided, since they
# concentrate the extreme samples in one of the two sets.
cat(">>> Splitting data (70% train / 30% test, stratified random)...\n")
set.seed(1234)
idx     <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[idx, ];  X_test  <- X[-idx, ]
y_train <- y[idx];    y_test  <- y[-idx]
cat("   Train:", nrow(X_train), "| Test:", nrow(X_test), "\n")
cat("   Train class dist:", paste(names(table(y_train)), table(y_train), sep = "=", collapse = " / "), "\n")
cat("   Test  class dist:", paste(names(table(y_test)),  table(y_test),  sep = "=", collapse = " / "), "\n\n")

# ======================================
# 6. Helper functions
# ======================================
apply_preprocessing <- function(Xtr, Xte, ytr, yte, scatter, deriv) {
  tryCatch({
    Xtr <- as.matrix(Xtr); Xte <- as.matrix(Xte)
    co  <- colnames(Xtr); wo <- parse_wn(co)

    # --- Order of operations (SCATTER_FIRST / SMOOTH_ALWAYS switches) ---
    # Everything is computed on the FULL, continuous spectrum (545 variables, no
    # gaps): this way the 11-point Savitzky-Golay moving window does not mix signal
    # from non-contiguous regions. Trimming to the Patz windows is done afterwards.
    keep_names <- function(Xnew, cn) {
      # Savitzky-Golay may trim (w-1)/2 columns at each edge; the remaining columns
      # are the central ones, so they are named by position.
      if (is.null(cn) || ncol(Xnew) == 0) return(Xnew)
      d <- length(cn) - ncol(Xnew)
      if (d >= 0 && d %% 2 == 0) colnames(Xnew) <- cn[(d / 2 + 1):(d / 2 + ncol(Xnew))]
      Xnew
    }
    sg_step <- function(Xa, Xb) {
      if (deriv > 0 || SMOOTH_ALWAYS) {
        cn <- colnames(Xa)
        Xa <- savitzkyGolay(Xa, m = deriv, p = 3, w = 11)
        Xb <- savitzkyGolay(Xb, m = deriv, p = 3, w = 11)
        Xa <- keep_names(Xa, cn); Xb <- keep_names(Xb, cn)
      }
      list(Xa, Xb)
    }
    scatter_step <- function(Xa, Xb) {
      cn <- colnames(Xa)
      if (scatter == "SNV") { Xa <- standardNormalVariate(Xa); Xb <- standardNormalVariate(Xb) }
      if (scatter == "MSC") {
        ref <- colMeans(Xa, na.rm = TRUE)   # reference: training set mean
        mf  <- function(r) { if (all(is.na(r))) return(r); fit <- lm(r ~ ref); (r - coef(fit)[1]) / coef(fit)[2] }
        Xa <- t(apply(Xa, 1, mf)); Xb <- t(apply(Xb, 1, mf))
      }
      if (!is.null(cn) && ncol(Xa) == length(cn)) { colnames(Xa) <- cn; colnames(Xb) <- cn }
      list(Xa, Xb)
    }
    if (SCATTER_FIRST) {
      s1 <- scatter_step(Xtr, Xte); s2 <- sg_step(s1[[1]], s1[[2]])
    } else {
      s1 <- sg_step(Xtr, Xte);      s2 <- scatter_step(s1[[1]], s1[[2]])
    }
    Xtr <- s2[[1]]; Xte <- s2[[2]]

    # --- Trim to the Patz windows AFTER smoothing/derivatization/scatter correction ---
    co2 <- colnames(Xtr); wo2 <- suppressWarnings(parse_wn(co2))
    if (length(wo2) == 0 || all(is.na(wo2))) wo2 <- wo
    keep_patz <- in_patz_windows(wo2)
    Xtr <- Xtr[, keep_patz, drop = FALSE]; Xte <- Xte[, keep_patz, drop = FALSE]
    wo2 <- wo2[keep_patz]

    mu  <- colMeans(Xtr, na.rm = TRUE)
    Xtr_plot <- Xtr; Xte_plot <- Xte   # <- uncentered version, for diagnostic plots only
    Xtr <- sweep(Xtr, 2, mu, "-"); Xte <- sweep(Xte, 2, mu, "-")
    ktr <- complete.cases(Xtr); kte <- complete.cases(Xte)
    vv  <- apply(Xtr[ktr, , drop = FALSE], 2, var, na.rm = TRUE)
    kc  <- !is.na(vv) & vv > 1e-10
    fn  <- colnames(Xtr)[kc]; fw <- suppressWarnings(parse_wn(fn))
    if (any(is.na(fw))) fw <- wo2[kc]
    list(Xtr = Xtr[ktr, kc, drop = FALSE], Xte = Xte[kte, kc, drop = FALSE],
         Xtr_plot = Xtr_plot[ktr, kc, drop = FALSE],
         ytr = ytr[ktr], yte = yte[kte], wavelengths = fw, success = TRUE)
  }, error = function(e) list(success = FALSE, error = as.character(e)))
}

plot_spectra <- function(Xm, wl, title, sel_vars = NULL, pp_name) {
  nc <- ncol(Xm)
  if (length(wl) > nc) wl <- wl[1:nc] else if (length(wl) < nc) wl <- 1:nc
  if (nrow(Xm) == 0 || nc == 0) return(NULL)
  df <- data.frame(wavelength = rep(wl, each = nrow(Xm)),
                    absorbance = as.vector(Xm),
                    sample     = rep(1:nrow(Xm), times = length(wl)))
  df <- df[complete.cases(df), ]; if (nrow(df) == 0) return(NULL)
  p  <- ggplot(df, aes(x = wavelength, y = absorbance, group = sample)) +
    geom_line(alpha = 0.3, color = "gray40") +
    labs(title = title, x = expression("Wavenumber (cm"^{-1}*")"), y = "Absorbance") +
    theme_minimal(base_size = 15) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          panel.grid.minor = element_blank())
  if (!is.null(sel_vars)) {
    sw <- parse_wn(sel_vars); sw <- sw[sw %in% wl]
    if (length(sw) > 0)
      p <- p + geom_vline(xintercept = sw, color = "red", alpha = 0.6, linetype = "dashed") +
        annotate("text", x = min(wl), y = max(df$absorbance, na.rm = TRUE),
                 label = paste0(length(sw), " variables selected"),
                 hjust = 0, vjust = 1, color = "red", size = 5, fontface = "bold")
  }
  fn <- paste0("Results/Spectra/", gsub(" ", "_", pp_name),
               ifelse(is.null(sel_vars), "", "_Boruta"), ".png")
  tryCatch(ggsave(fn, p, width = 10, height = 6, dpi = 300, bg = "white"), error = function(e) NULL)
  invisible(p)
}

# Complete binary metrics
calc_class_metrics <- function(model, tr, te, is_xgb = FALSE,
                                xgb_model = NULL, xgb_tr = NULL, xgb_te = NULL) {
  tryCatch({
    if (is_xgb) {
      # Predictions for manual XGB
      pred_tr_prob <- predict(xgb_model, xgb_tr)
      pred_te_prob <- predict(xgb_model, xgb_te)
      lvls         <- class_names
      pred_tr      <- factor(ifelse(pred_tr_prob > 0.5, lvls[2], lvls[1]), levels = lvls)
      pred_te      <- factor(ifelse(pred_te_prob > 0.5, lvls[2], lvls[1]), levels = lvls)
      y_tr         <- tr$Clase
      y_te         <- te$Clase
    } else {
      pred_tr <- predict(model, tr)
      pred_te <- predict(model, te)
      y_tr    <- tr$Clase
      y_te    <- te$Clase
    }

    cm_tr <- confusionMatrix(pred_tr, y_tr, positive = POS_CLASS)
    cm_te <- confusionMatrix(pred_te, y_te, positive = POS_CLASS)

    s_tr  <- cm_tr$byClass["Sensitivity"]; p_tr <- cm_tr$byClass["Pos Pred Value"]
    f1_tr <- ifelse(is.na(s_tr + p_tr) || (s_tr + p_tr) == 0, NA, 2 * s_tr * p_tr / (s_tr + p_tr))
    s_te  <- cm_te$byClass["Sensitivity"]; p_te <- cm_te$byClass["Pos Pred Value"]
    f1_te <- ifelse(is.na(s_te + p_te) || (s_te + p_te) == 0, NA, 2 * s_te * p_te / (s_te + p_te))

    auc_tr <- NA_real_; auc_te <- NA_real_
    tryCatch({
      if (is_xgb) {
        prob_tr <- pred_tr_prob
        prob_te <- pred_te_prob
      } else {
        prob_tr <- predict(model, tr, type = "prob")[, POS_CLASS]
        prob_te <- predict(model, te, type = "prob")[, POS_CLASS]
      }
      auc_tr <- as.numeric(roc(y_tr, prob_tr, levels = class_names, direction = "<", quiet = TRUE)$auc)
      auc_te <- as.numeric(roc(y_te, prob_te, levels = class_names, direction = "<", quiet = TRUE)$auc)
    }, error = function(e) NULL)

    list(Accuracy_Train    = cm_tr$overall["Accuracy"],    Accuracy_Test    = cm_te$overall["Accuracy"],
         Kappa_Train       = cm_tr$overall["Kappa"],        Kappa_Test       = cm_te$overall["Kappa"],
         Sensitivity_Train = cm_tr$byClass["Sensitivity"],  Sensitivity_Test = cm_te$byClass["Sensitivity"],
         Specificity_Train = cm_tr$byClass["Specificity"],  Specificity_Test = cm_te$byClass["Specificity"],
         Precision_Train   = cm_tr$byClass["Pos Pred Value"], Precision_Test = cm_te$byClass["Pos Pred Value"],
         F1_Train = f1_tr, F1_Test = f1_te,
         AUC_Train = auc_tr, AUC_Test = auc_te,
         CM_Train = cm_tr, CM_Test = cm_te, success = TRUE)
  }, error = function(e) list(success = FALSE, error = as.character(e)))
}

# Readable axis labels (dots -> spaces, line break)
wrap_cls <- function(x) {
  vapply(strwrap(gsub(".", " ", x, fixed = TRUE), width = 14, simplify = FALSE),
         paste, character(1), collapse = "\n")
}

# Confusion matrix plot
plot_confusion <- function(cm, model_name, pp_name, boruta_status, set_name) {
  tryCatch({
    tbl <- as.data.frame(cm$table)
    colnames(tbl) <- c("Prediction", "Reference", "Freq")   # as.data.frame(cm$table): 1st dim = Prediction, 2nd = Reference
    acc <- round(cm$overall["Accuracy"] * 100, 1); kap <- round(cm$overall["Kappa"], 3)
    p   <- ggplot(tbl, aes(x = Reference, y = Prediction, fill = Freq)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = Freq), color = "black", fontface = "bold", size = 7) +
      scale_fill_gradientn(colours = c("#f7fbff", "#2171b5", "#08306b"), name = "Count") +
      scale_x_discrete(labels = wrap_cls) + scale_y_discrete(labels = wrap_cls) +
      labs(title    = paste0("Confusion Matrix — ", model_name, " (", set_name, ")"),
           subtitle = paste0(pp_name, " | Boruta: ", boruta_status,
                             "  |  Acc=", acc, "%  |  κ=", kap),
           x = "Actual class", y = "Predicted class") +
      theme_minimal(base_size = 16) +
      theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 17),
            plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
            panel.grid    = element_blank(),
            axis.title    = element_text(size = 15),
            axis.text     = element_text(size = 13, face = "bold"))
    dir.create(paste0("Results/ConfusionMatrices/", model_name), showWarnings = FALSE, recursive = TRUE)
    fn <- paste0("Results/ConfusionMatrices/", model_name, "/",
                 gsub(" ", "_", pp_name), "_", boruta_status, "_", set_name, ".png")
    ggsave(fn, p, width = 7, height = 5.5, dpi = 300, bg = "white")
  }, error = function(e) NULL)
}

# Binary ROC curve
plot_roc_binary <- function(prob_vec, y_vec, model_name, pp_name, boruta_status, auc_val) {
  tryCatch({
    roc_obj <- roc(y_vec, prob_vec, levels = class_names, direction = "<", quiet = TRUE)
    roc_df  <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
    p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_ribbon(aes(ymin = 0, ymax = TPR), fill = "#377EB8", alpha = 0.15) +
      geom_line(color = "#377EB8", linewidth = 1.4) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
      annotate("text", x = 0.65, y = 0.12, label = paste0("AUC = ", round(auc_val, 4)),
               size = 6, fontface = "bold", color = "#377EB8") +
      scale_x_continuous(labels = percent, limits = c(0, 1)) +
      scale_y_continuous(labels = percent, limits = c(0, 1)) +
      labs(title    = paste0("ROC Curve — ", model_name),
           subtitle = paste0(pp_name, " | Boruta: ", boruta_status, " | Positive class: ", POS_CLASS),
           x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)") +
      theme_minimal(base_size = 16) +
      theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 17),
            plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
            axis.title    = element_text(size = 15), axis.text = element_text(size = 13),
            panel.border  = element_rect(color = "gray70", fill = NA))
    dir.create(paste0("Results/ROC/", model_name), showWarnings = FALSE, recursive = TRUE)
    fn <- paste0("Results/ROC/", model_name, "/",
                 gsub(" ", "_", pp_name), "_", boruta_status, ".png")
    ggsave(fn, p, width = 7, height = 6, dpi = 300, bg = "white")
  }, error = function(e) NULL)
}

# ======================================================
# HEATMAP TYPE 1: facet_grid
# ======================================================
make_heatmap_facet <- function(df_wide, x_var, y_var = "Model",
                                title, subtitle, filename,
                                width = 18, height = 11) {
  lbl_acc <- "Accuracy Test"
  lbl_f1  <- "F1 Test"
  lbl_auc <- "AUC Test"

  df_long <- df_wide %>%
    rename(x_col = !!sym(x_var), Model = !!sym(y_var)) %>%
    pivot_longer(cols = c(Accuracy_Test, F1_Test, AUC_Test),
                 names_to = "metric_id", values_to = "Value") %>%
    mutate(Metric = factor(
      case_when(metric_id == "Accuracy_Test" ~ lbl_acc,
                metric_id == "F1_Test"       ~ lbl_f1,
                metric_id == "AUC_Test"      ~ lbl_auc),
      levels = c(lbl_acc, lbl_f1, lbl_auc)))

  d_acc <- filter(df_long, metric_id == "Accuracy_Test")
  d_f1  <- filter(df_long, metric_id == "F1_Test")
  d_auc <- filter(df_long, metric_id == "AUC_Test")

  p <- ggplot() +
    geom_tile(data = d_acc, aes(x = x_col, y = Metric, fill = Value),
              color = "white", linewidth = 0.5) +
    geom_text(data = d_acc, aes(x = x_col, y = Metric, label = paste0(round(Value * 100, 1), "%")),
              color = "black", fontface = "bold", size = 2.7) +
    scale_fill_gradientn(
      name    = "Accuracy",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.75, 0.90, 1.0)), limits = c(0, 1),
      labels  = percent,
      guide   = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.4, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 8, face = "bold", color = "black"),
                               label.theme = element_text(size = 7, color = "black"), order = 1)) +
    new_scale_fill() +
    geom_tile(data = d_f1, aes(x = x_col, y = Metric, fill = Value),
              color = "white", linewidth = 0.5) +
    geom_text(data = d_f1, aes(x = x_col, y = Metric, label = round(Value, 3)),
              color = "black", fontface = "bold", size = 2.7) +
    scale_fill_gradientn(
      name    = "F1 Score",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.75, 0.90, 1.0)), limits = c(0, 1),
      guide   = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.4, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 8, face = "bold", color = "black"),
                               label.theme = element_text(size = 7, color = "black"), order = 2)) +
    new_scale_fill() +
    geom_tile(data = d_auc, aes(x = x_col, y = Metric, fill = Value),
              color = "white", linewidth = 0.5) +
    geom_text(data = d_auc, aes(x = x_col, y = Metric,
              label = ifelse(is.na(Value), "N/A", round(Value, 3))),
              color = "black", fontface = "bold", size = 2.7) +
    scale_fill_gradientn(
      name    = "AUC",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.80, 0.95, 1.0)), limits = c(0, 1),
      guide   = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.4, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 8, face = "bold", color = "black"),
                               label.theme = element_text(size = 7, color = "black"), order = 3)) +
    facet_grid(rows = vars(Model), scales = "free_y", space = "free_y", switch = "y") +
    labs(title = title, subtitle = subtitle, x = x_var, y = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 9.5) +
    theme(
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle     = element_text(hjust = 0.5, size = 7.5, color = "gray40", margin = margin(b = 6)),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9,
                                        hjust = 0.5, vjust = 0.5, margin = margin(r = 5, l = 5)),
      strip.background  = element_rect(fill = "gray88", color = "black", linewidth = 0.9),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.9),
      panel.spacing.y   = unit(2, "pt"),
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 7.5, color = "black"),
      axis.text.y       = element_text(size = 8, color = "black"),
      axis.ticks        = element_blank(),
      panel.grid        = element_blank(),
      plot.margin       = margin(8, 5, 8, 5),
      legend.position   = "right",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.5, "cm"),
      legend.key.size   = unit(0.4, "cm"))

  ggsave(filename, p, width = width, height = height, dpi = 300, bg = "white")
  cat("   OK Heatmap saved:", basename(filename), "\n")
  invisible(p)
}

# ======================================================
# HEATMAP TYPE 2: simple table
# ======================================================
make_heatmap_simple <- function(df_wide, row_var,
                                 title, subtitle, filename,
                                 width = 10, height = 7) {
  lbl_acc <- "Accuracy\nTest"
  lbl_f1  <- "F1\nTest"
  lbl_auc <- "AUC\nTest"

  df_long <- df_wide %>%
    rename(Row = !!sym(row_var)) %>%
    pivot_longer(cols = c(Accuracy_Test, F1_Test, AUC_Test),
                 names_to = "metric_id", values_to = "Value") %>%
    mutate(Metric = factor(
      case_when(metric_id == "Accuracy_Test" ~ lbl_acc,
                metric_id == "F1_Test"       ~ lbl_f1,
                metric_id == "AUC_Test"      ~ lbl_auc),
      levels = c(lbl_acc, lbl_f1, lbl_auc)))

  d_acc <- filter(df_long, metric_id == "Accuracy_Test")
  d_f1  <- filter(df_long, metric_id == "F1_Test")
  d_auc <- filter(df_long, metric_id == "AUC_Test")

  p <- ggplot() +
    geom_tile(data = d_acc, aes(x = Metric, y = Row, fill = Value), color = "white", linewidth = 0.7) +
    geom_text(data = d_acc, aes(x = Metric, y = Row, label = paste0(round(Value * 100, 1), "%")),
              color = "black", fontface = "bold", size = 3.5) +
    scale_fill_gradientn(name = "Accuracy",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.75, 0.90, 1.0)), limits = c(0, 1), labels = percent,
      guide   = guide_colorbar(barheight = unit(3.5, "cm"), barwidth = unit(0.5, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 9, face = "bold", color = "black"),
                               label.theme = element_text(size = 8, color = "black"), order = 1)) +
    new_scale_fill() +
    geom_tile(data = d_f1, aes(x = Metric, y = Row, fill = Value), color = "white", linewidth = 0.7) +
    geom_text(data = d_f1, aes(x = Metric, y = Row, label = round(Value, 3)),
              color = "black", fontface = "bold", size = 3.5) +
    scale_fill_gradientn(name = "F1 Score",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.75, 0.90, 1.0)), limits = c(0, 1),
      guide   = guide_colorbar(barheight = unit(3.5, "cm"), barwidth = unit(0.5, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 9, face = "bold", color = "black"),
                               label.theme = element_text(size = 8, color = "black"), order = 2)) +
    new_scale_fill() +
    geom_tile(data = d_auc, aes(x = Metric, y = Row, fill = Value), color = "white", linewidth = 0.7) +
    geom_text(data = d_auc, aes(x = Metric, y = Row,
              label = ifelse(is.na(Value), "N/A", round(Value, 3))),
              color = "black", fontface = "bold", size = 3.5) +
    scale_fill_gradientn(name = "AUC",
      colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
      values  = rescale(c(0, 0.60, 0.80, 0.95, 1.0)), limits = c(0, 1),
      guide   = guide_colorbar(barheight = unit(3.5, "cm"), barwidth = unit(0.5, "cm"),
                               title.position = "top", title.hjust = 0.5,
                               title.theme = element_text(size = 9, face = "bold", color = "black"),
                               label.theme = element_text(size = 8, color = "black"), order = 3)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40", margin = margin(b = 6)),
          axis.text.x   = element_text(angle = 0, hjust = 0.5, size = 10, face = "bold", color = "black"),
          axis.text.y   = element_text(size = 10, color = "black"),
          axis.ticks    = element_blank(), panel.grid = element_blank(),
          legend.position  = "right", legend.box = "vertical",
          legend.spacing.y = unit(0.5, "cm"), legend.key.size = unit(0.5, "cm"),
          plot.margin   = margin(8, 5, 8, 5))

  ggsave(filename, p, width = width, height = height, dpi = 300, bg = "white")
  cat("   OK Heatmap saved:", basename(filename), "\n")
  invisible(p)
}

# ======================================
# 7. Classification models
# ======================================
ctrl_cls <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  allowParallel   = TRUE,
  verboseIter     = FALSE,
  savePredictions = "final"
)

# ── XGB: helper function that trains directly with xgboost ────────
# Completely avoids the caret wrapper that causes "Error: Stopping"
train_xgb_direct <- function(tr) {
  set.seed(123)

  lvls    <- levels(tr$Clase)
  y_num   <- as.integer(tr$Clase) - 1   # 0/1
  X_mat   <- as.matrix(tr[, setdiff(names(tr), "Clase")])
  dtrain  <- xgb.DMatrix(data = X_mat, label = y_num)

  # Manual search grid (36 reasonable combinations)
  grid <- expand.grid(
    nrounds          = c(50, 100, 150),
    max_depth        = c(2, 3, 4),
    eta              = c(0.05, 0.1),
    gamma            = c(0, 0.1),
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample        = 0.8,
    stringsAsFactors = FALSE
  )

  best_auc   <- -Inf
  best_model <- NULL
  best_params <- NULL

  nfolds <- 5
  folds  <- createFolds(tr$Clase, k = nfolds, list = TRUE)

  for (i in seq_len(nrow(grid))) {
    params <- list(
      objective        = "binary:logistic",
      eval_metric      = "auc",
      max_depth        = grid$max_depth[i],
      eta              = grid$eta[i],
      gamma            = grid$gamma[i],
      colsample_bytree = grid$colsample_bytree[i],
      min_child_weight = grid$min_child_weight[i],
      subsample        = grid$subsample[i],
      nthread          = 1,
      verbosity        = 0
    )

    # Manual CV
    cv_aucs <- numeric(nfolds)
    for (f in seq_len(nfolds)) {
      val_idx   <- folds[[f]]
      tr_idx    <- setdiff(seq_len(nrow(X_mat)), val_idx)
      d_tr      <- xgb.DMatrix(data = X_mat[tr_idx, , drop = FALSE], label = y_num[tr_idx])
      d_val     <- xgb.DMatrix(data = X_mat[val_idx, , drop = FALSE], label = y_num[val_idx])
      m_tmp     <- xgb.train(params = params, data = d_tr,
                              nrounds = grid$nrounds[i], verbose = 0)
      prob_val  <- predict(m_tmp, d_val)
      tryCatch({
        cv_aucs[f] <- as.numeric(roc(y_num[val_idx], prob_val, quiet = TRUE)$auc)
      }, error = function(e) { cv_aucs[f] <<- 0.5 })
    }

    mean_auc <- mean(cv_aucs, na.rm = TRUE)
    if (mean_auc > best_auc) {
      best_auc    <- mean_auc
      best_params <- list(params = params, nrounds = grid$nrounds[i])
    }
  }

  # Train the final model on all the data
  final_model <- xgb.train(
    params  = best_params$params,
    data    = dtrain,
    nrounds = best_params$nrounds,
    verbose = 0
  )

  list(model = final_model, levels = lvls, cv_auc = best_auc, best_params = best_params)
}

# Prediction for direct XGB
predict_xgb <- function(xgb_obj, newdata_df) {
  X_mat <- as.matrix(newdata_df[, setdiff(names(newdata_df), "Clase")])
  dmat  <- xgb.DMatrix(data = X_mat)
  probs <- predict(xgb_obj$model, dmat)
  list(prob = probs, class = factor(
    ifelse(probs > 0.5, xgb_obj$levels[2], xgb_obj$levels[1]),
    levels = xgb_obj$levels))
}

models <- list(

  LDA = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "lda", trControl = ctrl_cls,
          metric = "ROC", preProcess = c("center", "scale"))
  },

  PLSDA = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "pls", trControl = ctrl_cls,
          tuneLength = 15, metric = "ROC", preProcess = c("center", "scale"))
  },

  SVM = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "svmRadial", trControl = ctrl_cls,
          preProcess = c("center", "scale"), tuneLength = 8, metric = "ROC")
  },

  RF = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "rf", trControl = ctrl_cls,
          tuneLength = 6, ntree = 300, metric = "ROC")
  },

  DT = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "rpart", trControl = ctrl_cls,
          tuneGrid = data.frame(cp = 10^seq(-4, -1, length = 20)), metric = "ROC")
  },

  kNN = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "knn", trControl = ctrl_cls,
          preProcess = c("center", "scale"),
          tuneGrid   = data.frame(k = c(3, 5, 7, 9, 11, 13, 15)), metric = "ROC")
  },

  NB = function(tr) {
    set.seed(123)
    train(Clase ~ ., data = tr, method = "naive_bayes", trControl = ctrl_cls,
          tuneGrid = expand.grid(laplace = c(0, 1), usekernel = TRUE, adjust = c(0.5, 1, 1.5)),
          metric = "ROC")
  }
)

# ======================================
# 8. Experimental setup
# ======================================
scatter_opts <- c("none", "SNV", "MSC")
deriv_opts   <- c(0, 1, 2)
boruta_opts  <- c("none", "boruta")

results_all <- list(); error_log <- list(); boruta_vars_log <- list()
start_time <- Sys.time()

cat("\n>>> EXPERIMENT CONFIGURATION:\n")
cat("   Preprocessing methods:", length(scatter_opts) * length(deriv_opts), "\n")
cat("   Models:", length(models) + 1, "(including direct XGB)\n")
cat("   Boruta options:", length(boruta_opts), "\n")
cat("   TOTAL COMBINATIONS:",
    length(scatter_opts) * length(deriv_opts) * length(boruta_opts) * (length(models) + 1), "\n")
cat("   Estimated time: ~2-5 hours\n\n")

# ======================================
# 9. Main loop
# ======================================
total_runs  <- length(scatter_opts) * length(deriv_opts) * length(boruta_opts) * (length(models) + 1)
current_run <- 0

for (sc in scatter_opts) {
  for (dv in deriv_opts) {
    pp_name <- make_pp_name(sc, dv)
    cat("\n>>> Preprocessing:", pp_name, "\n")

    pp <- apply_preprocessing(X_train, X_test, y_train, y_test, sc, dv)
    if (!pp$success) {
      cat("   ERROR:", pp$error, "\n")
      error_log[[length(error_log) + 1]] <- list(step = pp_name, error = pp$error)
      next
    }

    plot_spectra(pp$Xtr_plot, pp$wavelengths, paste("Preprocessed spectra (uncentered):", pp_name), NULL, pp_name)

    for (sel in boruta_opts) {
      train_base <- data.frame(Clase = pp$ytr, pp$Xtr)
      test_base  <- data.frame(Clase = pp$yte, pp$Xte)

      if (sel == "boruta") {
        cat("   Running Boruta...\n")
        bor_res <- tryCatch({
          bor <- Boruta(Clase ~ ., train_base, maxRuns = 500, doTrace = 0)
          n_tentative <- sum(bor$finalDecision == "Tentative")
          if (n_tentative > 0) {
            cat("   ", n_tentative, "variable(s) Tentative -> resolving with TentativeRoughFix...\n")
            bor <- TentativeRoughFix(bor)
          }
          vars <- names(bor$finalDecision[bor$finalDecision == "Confirmed"])
          if (length(vars) == 0) {
            cat("   WARNING: no vars confirmed. Using all.\n")
            list(success = TRUE, use_all = TRUE, n_tentative = n_tentative)
          } else {
            cat("   Boruta selected", length(vars), "vars\n")
            list(success = TRUE, use_all = FALSE, vars = vars, n_tentative = n_tentative)
          }
        }, error = function(e) list(success = FALSE, error = as.character(e)))

        if (!bor_res$success) {
          error_log[[length(error_log) + 1]] <- list(
            step = paste(pp_name, "Boruta"), error = bor_res$error)
          next
        }
        boruta_vars_log[[paste(pp_name, "boruta", sep = "__")]] <- list(
          Preprocessing = pp_name,
          N_selected    = if (bor_res$use_all) 0 else length(bor_res$vars),
          N_tentative   = bor_res$n_tentative,
          Variables     = if (bor_res$use_all) "" else paste(bor_res$vars, collapse = ", "),
          use_all       = bor_res$use_all)
        if (!bor_res$use_all) {
          train_base <- train_base[, c("Clase", bor_res$vars)]
          test_base  <- test_base[,  c("Clase", bor_res$vars)]
          plot_spectra(pp$Xtr_plot, pp$wavelengths, paste("Boruta vars:", pp_name), bor_res$vars, pp_name)
        }
      }

      # ── caret models ──────────────────────────────────────────────
      for (m in names(models)) {
        current_run <- current_run + 1
        cat("   [", current_run, "/", total_runs, "] Model:", m, "\n")

        res <- tryCatch({
          model <- models[[m]](train_base)
          ev    <- calc_class_metrics(model, train_base, test_base)
          if (!ev$success) return(list(success = FALSE, error = ev$error))

          plot_confusion(ev$CM_Train, m, pp_name, sel, "Train")
          plot_confusion(ev$CM_Test,  m, pp_name, sel, "Test")

          if (!is.na(ev$AUC_Test)) {
            prob_te <- predict(model, test_base, type = "prob")[, POS_CLASS]
            plot_roc_binary(prob_te, test_base$Clase, m, pp_name, sel, ev$AUC_Test)
          }
          list(success = TRUE, metrics = ev)
        }, error = function(e) list(success = FALSE, error = as.character(e)))

        if (res$success) {
          ev <- res$metrics
          results_all[[length(results_all) + 1]] <- data.frame(
            Preprocessing    = pp_name, Boruta = sel, Model = m, N_vars = ncol(train_base) - 1,
            Accuracy_Train   = ev$Accuracy_Train,    Accuracy_Test    = ev$Accuracy_Test,
            Kappa_Train      = ev$Kappa_Train,        Kappa_Test       = ev$Kappa_Test,
            Sensitivity_Train = ev$Sensitivity_Train, Sensitivity_Test = ev$Sensitivity_Test,
            Specificity_Train = ev$Specificity_Train, Specificity_Test = ev$Specificity_Test,
            Precision_Train  = ev$Precision_Train,    Precision_Test   = ev$Precision_Test,
            F1_Train  = ev$F1_Train,  F1_Test  = ev$F1_Test,
            AUC_Train = ev$AUC_Train, AUC_Test = ev$AUC_Test,
            Hyperparameters = format_bestTune(model$bestTune))
        } else {
          cat("   ERROR:", res$error, "\n")
          error_log[[length(error_log) + 1]] <- list(step = paste(pp_name, sel, m), error = res$error)
        }
      }

      # ── Direct XGB (without caret) ─────────────────────────────────────
      current_run <- current_run + 1
      cat("   [", current_run, "/", total_runs, "] Model: XGB\n")

      res_xgb <- tryCatch({
        xgb_obj  <- train_xgb_direct(train_base)
        pred_tr  <- predict_xgb(xgb_obj, train_base)
        pred_te  <- predict_xgb(xgb_obj, test_base)

        y_tr_f <- factor(train_base$Clase, levels = class_names)
        y_te_f <- factor(test_base$Clase,  levels = class_names)

        cm_tr <- confusionMatrix(pred_tr$class, y_tr_f, positive = POS_CLASS)
        cm_te <- confusionMatrix(pred_te$class, y_te_f, positive = POS_CLASS)

        s_tr  <- cm_tr$byClass["Sensitivity"]; p_tr_v <- cm_tr$byClass["Pos Pred Value"]
        f1_tr <- ifelse(is.na(s_tr + p_tr_v) || (s_tr + p_tr_v) == 0, NA,
                        2 * s_tr * p_tr_v / (s_tr + p_tr_v))
        s_te  <- cm_te$byClass["Sensitivity"]; p_te_v <- cm_te$byClass["Pos Pred Value"]
        f1_te <- ifelse(is.na(s_te + p_te_v) || (s_te + p_te_v) == 0, NA,
                        2 * s_te * p_te_v / (s_te + p_te_v))

        auc_tr <- tryCatch(
          as.numeric(roc(as.integer(y_tr_f) - 1, pred_tr$prob, quiet = TRUE)$auc),
          error = function(e) NA_real_)
        auc_te <- tryCatch(
          as.numeric(roc(as.integer(y_te_f) - 1, pred_te$prob, quiet = TRUE)$auc),
          error = function(e) NA_real_)

        plot_confusion(cm_tr, "XGB", pp_name, sel, "Train")
        plot_confusion(cm_te, "XGB", pp_name, sel, "Test")
        if (!is.na(auc_te))
          plot_roc_binary(pred_te$prob, y_te_f, "XGB", pp_name, sel, auc_te)

        list(success = TRUE,
             Accuracy_Train    = cm_tr$overall["Accuracy"],    Accuracy_Test    = cm_te$overall["Accuracy"],
             Kappa_Train       = cm_tr$overall["Kappa"],        Kappa_Test       = cm_te$overall["Kappa"],
             Sensitivity_Train = cm_tr$byClass["Sensitivity"],  Sensitivity_Test = cm_te$byClass["Sensitivity"],
             Specificity_Train = cm_tr$byClass["Specificity"],  Specificity_Test = cm_te$byClass["Specificity"],
             Precision_Train   = cm_tr$byClass["Pos Pred Value"], Precision_Test = cm_te$byClass["Pos Pred Value"],
             F1_Train = f1_tr, F1_Test = f1_te,
             AUC_Train = auc_tr, AUC_Test = auc_te,
             Hyperparameters = format_xgb_params(xgb_obj$best_params))
      }, error = function(e) list(success = FALSE, error = as.character(e)))

      if (res_xgb$success) {
        results_all[[length(results_all) + 1]] <- data.frame(
          Preprocessing     = pp_name, Boruta = sel, Model = "XGB", N_vars = ncol(train_base) - 1,
          Accuracy_Train    = res_xgb$Accuracy_Train,    Accuracy_Test    = res_xgb$Accuracy_Test,
          Kappa_Train       = res_xgb$Kappa_Train,        Kappa_Test       = res_xgb$Kappa_Test,
          Sensitivity_Train = res_xgb$Sensitivity_Train,  Sensitivity_Test = res_xgb$Sensitivity_Test,
          Specificity_Train = res_xgb$Specificity_Train,  Specificity_Test = res_xgb$Specificity_Test,
          Precision_Train   = res_xgb$Precision_Train,    Precision_Test   = res_xgb$Precision_Test,
          F1_Train  = res_xgb$F1_Train,  F1_Test  = res_xgb$F1_Test,
          AUC_Train = res_xgb$AUC_Train, AUC_Test = res_xgb$AUC_Test,
          Hyperparameters = res_xgb$Hyperparameters)
      } else {
        cat("   ERROR XGB:", res_xgb$error, "\n")
        error_log[[length(error_log) + 1]] <- list(step = paste(pp_name, sel, "XGB"), error = res_xgb$error)
      }
    }
  }
}

# ======================================
# 10. Export results and plots
# ======================================
if (length(results_all) > 0) {
  results_df <- bind_rows(results_all)
  write.xlsx(results_df, "Results/Excel/Complete_Summary.xlsx", overwrite = TRUE)
  cat("\nOK Results exported:", nrow(results_df), "rows\n")

  cat("\n=== TOP 10 (AUC Test) ===\n")
  print(
    results_df %>%
      arrange(desc(AUC_Test), desc(Accuracy_Test), desc(Kappa_Test), N_vars) %>%
      head(10) %>%
      dplyr::select(Model, Preprocessing, Boruta, N_vars, AUC_Test, Accuracy_Test, F1_Test, Kappa_Test)
  )

  # ────────────────────────────────────────────────────────────────
  # 9bis. Export the optimized hyperparameters of the top 10 models
  # (by AUC_Test) to Excel, together with their metrics.
  # ────────────────────────────────────────────────────────────────
  top10_hp <- results_df %>%
    arrange(desc(AUC_Test), desc(Accuracy_Test), desc(Kappa_Test), N_vars) %>%
    head(10) %>%
    mutate(Rank = row_number(), .before = 1) %>%
    dplyr::select(Rank, Model, Preprocessing, Boruta, N_vars, Hyperparameters,
                  AUC_Test, Accuracy_Test, F1_Test, Kappa_Test,
                  AUC_Train, Accuracy_Train, F1_Train, Kappa_Train)
  write.xlsx(top10_hp, "Results/Excel/Top10_Hyperparameters.xlsx", overwrite = TRUE)
  cat("   OK Saved: Top10_Hyperparameters.xlsx\n")

  # ────────────────────────────────────────────────────────────────
  # 10bis. Export the variables selected by Boruta to Excel
  #   (a) All Preprocessing x Boruta combinations, one row
  #       per combination, with the full list of selected wavenumbers
  #       -> to assess the overall effect of Boruta.
  #   (b) Only the combination of the WINNING MODEL (highest AUC_Test among
  #       the rows with Boruta=="boruta"), in long format (one wavenumber
  #       per row) to subsequently cross-reference it with the EDTA
  #       functional groups.
  # ────────────────────────────────────────────────────────────────
  if (length(boruta_vars_log) > 0) {
    cat("\n>>> Exporting Boruta variable selection to Excel...\n")

    bvl_all_df <- bind_rows(lapply(names(boruta_vars_log), function(k) {
      x <- boruta_vars_log[[k]]
      data.frame(
        Config_Key    = k,
        Preprocessing = x$Preprocessing,
        N_selected    = x$N_selected,
        N_tentative   = x$N_tentative,
        Used_all_vars = isTRUE(x$use_all),
        Variables     = x$Variables,
        stringsAsFactors = FALSE)
    }))
    write.xlsx(bvl_all_df,
               "Results/Excel/Boruta_Variables_All_Preprocessing.xlsx",
               overwrite = TRUE)
    cat("   OK Saved: Boruta_Variables_All_Preprocessing.xlsx (",
        nrow(bvl_all_df), "preprocessing combinations )\n")

    best_boruta_row <- results_df %>%
      filter(Boruta == "boruta") %>%
      arrange(desc(AUC_Test), desc(Accuracy_Test), desc(Kappa_Test), N_vars) %>%
      dplyr::slice(1)

    if (nrow(best_boruta_row) == 1) {
      best_key <- paste(best_boruta_row$Preprocessing, "boruta", sep = "__")
      best_log <- boruta_vars_log[[best_key]]

      if (!is.null(best_log) && nzchar(best_log$Variables)) {
        best_vars_wn <- suppressWarnings(parse_wn(
          trimws(strsplit(best_log$Variables, ",")[[1]])))
        best_vars_wn <- sort(best_vars_wn[!is.na(best_vars_wn)])

        best_vars_df <- data.frame(
          Winning_Model   = best_boruta_row$Model,
          Preprocessing   = best_boruta_row$Preprocessing,
          N_vars_selected = length(best_vars_wn),
          Wavenumber_cm1  = best_vars_wn,
          AUC_Test        = round(best_boruta_row$AUC_Test, 3),
          Accuracy_Test   = round(best_boruta_row$Accuracy_Test, 3),
          F1_Test         = round(best_boruta_row$F1_Test, 3))

        write.xlsx(best_vars_df,
                   "Results/Excel/Boruta_Variables_Best_Model.xlsx",
                   overwrite = TRUE)
        cat("   OK Saved: Boruta_Variables_Best_Model.xlsx  (Model:",
            best_boruta_row$Model, "| Preprocessing:", best_boruta_row$Preprocessing,
            "|", length(best_vars_wn), "wavenumbers )\n")
      } else {
        cat("   WARNING: winning combination used ALL variables (Boruta found no Confirmed vars) or log missing (",
            best_key, "). Skipping Boruta_Variables_Best_Model.xlsx\n")
      }
    } else {
      cat("   WARNING: no rows with Boruta=='boruta' in results_df. Skipping Boruta_Variables_Best_Model.xlsx\n")
    }
  }

  subtitle_hm <- paste0(
    "Accuracy/F1: red<0.75, yellow 0.75\u20130.90, green\u22650.90  |  ",
    "AUC: red<0.80, yellow 0.80\u20130.95, green\u22650.95")

  cat("\n>>> Generating heatmaps...\n")

  hm1 <- results_df %>% group_by(Model, Preprocessing) %>%
    summarise(Accuracy_Test = mean(Accuracy_Test, na.rm = TRUE),
              F1_Test       = mean(F1_Test,       na.rm = TRUE),
              AUC_Test      = mean(AUC_Test,      na.rm = TRUE), .groups = "drop")
  make_heatmap_facet(hm1, "Preprocessing", "Model",
    "Classification Metrics (Test) \u2014 Model vs Preprocessing",
    subtitle_hm, "Results/Heatmaps/Metrics_Model_Preprocessing.png", width = 14, height = 11)

  hm2 <- results_df %>% mutate(Config = paste(Preprocessing, Boruta, sep = "\n")) %>%
    group_by(Model, Config) %>%
    summarise(Accuracy_Test = mean(Accuracy_Test, na.rm = TRUE),
              F1_Test       = mean(F1_Test,       na.rm = TRUE),
              AUC_Test      = mean(AUC_Test,      na.rm = TRUE), .groups = "drop")
  make_heatmap_facet(hm2, "Config", "Model",
    "Classification Metrics (Test) \u2014 Model vs Preprocessing + Boruta",
    subtitle_hm, "Results/Heatmaps/Metrics_Model_Boruta.png", width = 22, height = 11)

  hm3 <- results_df %>% group_by(Model) %>%
    summarise(Accuracy_Test = mean(Accuracy_Test, na.rm = TRUE),
              F1_Test       = mean(F1_Test,       na.rm = TRUE),
              AUC_Test      = mean(AUC_Test,      na.rm = TRUE), .groups = "drop")
  make_heatmap_simple(hm3, "Model",
    "Classification Metrics (Test) \u2014 Average by Model",
    "Average across all preprocessing methods and Boruta options",
    "Results/Heatmaps/Metrics_by_Model.png", width = 9, height = 6)

  hm4 <- results_df %>% group_by(Preprocessing) %>%
    summarise(Accuracy_Test = mean(Accuracy_Test, na.rm = TRUE),
              F1_Test       = mean(F1_Test,       na.rm = TRUE),
              AUC_Test      = mean(AUC_Test,      na.rm = TRUE), .groups = "drop")
  make_heatmap_simple(hm4, "Preprocessing",
    "Classification Metrics (Test) \u2014 Average by Preprocessing",
    "Average across all models and Boruta options",
    "Results/Heatmaps/Metrics_by_Preprocessing.png", width = 9, height = 7)

  hm5 <- results_df %>% group_by(Model, Preprocessing) %>%
    summarise(AUC_Test = mean(AUC_Test, na.rm = TRUE), .groups = "drop")
  p_hm5 <- ggplot(hm5, aes(x = Preprocessing, y = Model, fill = AUC_Test)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = ifelse(is.na(AUC_Test), "N/A", round(AUC_Test, 3))),
              color = "black", fontface = "bold", size = 3.2) +
    scale_fill_gradientn(name = "AUC Test",
                         colours = c("red", "orange", "yellow", "yellowgreen", "darkgreen"),
                         values  = rescale(c(0, 0.60, 0.80, 0.95, 1.0)), limits = c(0, 1),
                         guide   = guide_colorbar(barheight = 8, barwidth = 1)) +
    labs(title    = "AUC Test \u2014 Model vs Preprocessing",
         subtitle = "Red<0.80 | Yellow 0.80\u20130.95 | Green\u22650.95",
         x = "Preprocessing", y = "Model") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
          axis.text.x   = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  ggsave("Results/Heatmaps/AUC_Model_Preprocessing.png", p_hm5, width = 12, height = 8, dpi = 300, bg = "white")
  cat("   OK AUC-only heatmap saved\n")

  preproc_pal <- setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#A65628", "#F781BF", "#999999", "#00CED1")[seq_len(length(unique(results_df$Preprocessing)))],
    unique(results_df$Preprocessing))

  p_comp <- ggplot(results_df, aes(x = Model, y = AUC_Test, color = Preprocessing)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95, ymax = Inf,   fill = "green",  alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.80, ymax = 0.95,  fill = "orange", alpha = 0.06) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.80,  fill = "red",    alpha = 0.04) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkgreen",  linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "darkorange", linewidth = 0.5, alpha = 0.7) +
    geom_jitter(aes(shape = Boruta), width = 0.25, size = 2.5, alpha = 0.75) +
    stat_summary(aes(group = Model), fun = median, geom = "crossbar",
                 width = 0.5, linewidth = 0.6, fatten = 2, show.legend = FALSE) +
    scale_color_manual(values = preproc_pal) +
    scale_shape_manual(values = c("none" = 16, "boruta" = 17),
                       labels = c("none" = "Without Boruta", "boruta" = "With Boruta")) +
    annotate("text", x = 0.6, y = 0.963, label = "AUC \u2265 0.95",
             color = "darkgreen", size = 3, hjust = 0, fontface = "italic") +
    annotate("text", x = 0.6, y = 0.813, label = "AUC 0.80\u20130.95",
             color = "darkorange", size = 3, hjust = 0, fontface = "italic") +
    labs(title    = "AUC Test by Model and Preprocessing",
         subtitle = "Points = individual combinations  |  Bar = median  |  Shapes = Boruta",
         x = "Model", y = "AUC (Test set)", color = "Preprocessing", shape = "Variable selection") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
          axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "right", panel.grid.major.x = element_blank())
  ggsave("Results/Other/Model_Comparison_AUC.png", p_comp, width = 13, height = 7, dpi = 300, bg = "white")

  p_comp_acc <- ggplot(results_df, aes(x = Model, y = Accuracy_Test, color = Preprocessing)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.90, ymax = Inf,  fill = "green",  alpha = 0.04) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.75, ymax = 0.90, fill = "orange", alpha = 0.06) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.75, fill = "red",    alpha = 0.04) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "darkgreen",  linewidth = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0.75, linetype = "dashed", color = "darkorange", linewidth = 0.5, alpha = 0.7) +
    geom_jitter(aes(shape = Boruta), width = 0.25, size = 2.5, alpha = 0.75) +
    stat_summary(aes(group = Model), fun = median, geom = "crossbar",
                 width = 0.5, linewidth = 0.6, fatten = 2, show.legend = FALSE) +
    scale_color_manual(values = preproc_pal) +
    scale_shape_manual(values = c("none" = 16, "boruta" = 17),
                       labels = c("none" = "Without Boruta", "boruta" = "With Boruta")) +
    scale_y_continuous(labels = percent) +
    labs(title    = "Accuracy Test by Model and Preprocessing",
         subtitle = "Points = individual combinations  |  Bar = median  |  Shapes = Boruta",
         x = "Model", y = "Accuracy (Test set)", color = "Preprocessing", shape = "Variable selection") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
          axis.text.x   = element_text(angle = 45, hjust = 1),
          legend.position = "right", panel.grid.major.x = element_blank())
  ggsave("Results/Other/Model_Comparison_Accuracy.png", p_comp_acc, width = 13, height = 7, dpi = 300, bg = "white")
  cat("   OK Model comparison plots saved\n")

  # PCA of metrics
  cat("\n>>> Generating PCA of metrics...\n")
  pca_cols <- c("Accuracy_Train", "Kappa_Train", "Sensitivity_Train", "Specificity_Train", "F1_Train", "AUC_Train",
                "Accuracy_Test",  "Kappa_Test",  "Sensitivity_Test",  "Specificity_Test",  "F1_Test",  "AUC_Test")
  pca_mat_raw <- results_df %>%
    dplyr::select(all_of(pca_cols)) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

  # Remove zero-variance columns before scaling
  col_var   <- apply(pca_mat_raw, 2, var, na.rm = TRUE)
  pca_mat   <- scale(pca_mat_raw[, col_var > 1e-10, drop = FALSE])

  pca_res   <- prcomp(pca_mat, center = FALSE, scale. = FALSE)
  pca_sc    <- as.data.frame(pca_res$x[, 1:2]); colnames(pca_sc) <- c("PC1", "PC2")
  pca_sc$Model         <- results_df$Model
  pca_sc$Preprocessing <- results_df$Preprocessing
  pca_sc$Boruta        <- results_df$Boruta
  pca_sc$AUC_Test      <- results_df$AUC_Test
  ve <- summary(pca_res)$importance[2, 1:2] * 100

  p_pca <- ggplot(pca_sc, aes(x = PC1, y = PC2, color = Model, shape = Boruta, size = AUC_Test)) +
    geom_point(alpha = 0.75) +
    scale_size_continuous(name = "AUC\n(Test)", range = c(1.5, 5)) +
    scale_shape_manual(values = c("none" = 16, "boruta" = 15),
                       labels = c("none" = "Without Boruta", "boruta" = "With Boruta")) +
    labs(title    = "PCA of Classification Metrics (Train + Test)",
         subtitle = "Point size proportional to AUC Test",
         x = paste0("PC1 (", round(ve[1], 1), "%)"), y = paste0("PC2 (", round(ve[2], 1), "%)"),
         color = "Model", shape = "Variable Selection") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "right", panel.border = element_rect(color = "gray50", fill = NA))
  ggsave("Results/Other/PCA_Metrics.png", p_pca, width = 14, height = 8, dpi = 300, bg = "white")

  ld <- as.data.frame(pca_res$rotation[, 1:2]); ld$Variable <- rownames(ld)
  write.xlsx(list("Scores"   = pca_sc, "Loadings" = ld,
                  "Variance_Explained" = data.frame(
                    Component    = paste0("PC", seq_along(pca_res$sdev)),
                    Var_Explained = summary(pca_res)$importance[2, ] * 100,
                    Cumulative    = summary(pca_res)$importance[3, ] * 100)),
             "Results/Excel/PCA_Analysis.xlsx", overwrite = TRUE)

  # Ranking
  cat("\n>>> Generating ranking...\n")
  ranking <- results_df %>%
    arrange(desc(AUC_Test), desc(F1_Test), desc(Accuracy_Test)) %>%
    mutate(Rank = row_number()) %>%
    dplyr::select(Rank, Model, Preprocessing, Boruta, N_vars,
                  AUC_Test, Accuracy_Test, Kappa_Test, Sensitivity_Test, Specificity_Test, F1_Test,
                  AUC_Train, Accuracy_Train, Kappa_Train, Sensitivity_Train, Specificity_Train, F1_Train)

  ms <- results_df %>% group_by(Model) %>%
    summarise(N = n(), AUC_Mean = mean(AUC_Test, na.rm = TRUE), AUC_Max = max(AUC_Test, na.rm = TRUE),
              Acc_Mean = mean(Accuracy_Test, na.rm = TRUE), F1_Mean = mean(F1_Test, na.rm = TRUE),
              Kappa_Mean = mean(Kappa_Test, na.rm = TRUE), .groups = "drop") %>% arrange(desc(AUC_Mean))
  ps <- results_df %>% group_by(Preprocessing) %>%
    summarise(N = n(), AUC_Mean = mean(AUC_Test, na.rm = TRUE),
              Acc_Mean = mean(Accuracy_Test, na.rm = TRUE), F1_Mean = mean(F1_Test, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(AUC_Mean))
  bs <- results_df %>% group_by(Boruta) %>%
    summarise(N = n(), AUC_Mean = mean(AUC_Test, na.rm = TRUE),
              Acc_Mean = mean(Accuracy_Test, na.rm = TRUE), F1_Mean = mean(F1_Test, na.rm = TRUE), .groups = "drop")

  # ======================================
  # 8bis. Paired t-test: effect of Boruta (AUC, Accuracy)
  #   Same approach as in the regression pipeline: one-tailed paired t-test,
  #   paired by Model x Preprocessing combination
  #   (same wine samples, different variable set).
  #   H1: AUC(Boruta) > AUC(None); Accuracy(Boruta) > Accuracy(None).
  # ======================================
  cat("\n>>> Paired t-test: Boruta effect (AUC, Accuracy)...\n")
  boruta_lv <- sort(unique(results_df$Boruta))
  bor_lab   <- boruta_lv[grepl("bor", boruta_lv, ignore.case = TRUE)]
  none_lab  <- boruta_lv[!grepl("bor", boruta_lv, ignore.case = TRUE)]

  if (length(boruta_lv) == 2 && length(bor_lab) == 1 && length(none_lab) == 1) {

    wide_metric <- function(metric_col) {
      results_df %>%
        select(Model, Preprocessing, Boruta, all_of(metric_col)) %>%
        tidyr::pivot_wider(names_from = Boruta, values_from = all_of(metric_col))
    }
    w_auc <- wide_metric("AUC_Test")
    w_acc <- wide_metric("Accuracy_Test")

    ok_auc <- complete.cases(w_auc[, c(bor_lab, none_lab)])
    ok_acc <- complete.cases(w_acc[, c(bor_lab, none_lab)])

    pt_auc <- t.test(w_auc[[bor_lab]][ok_auc], w_auc[[none_lab]][ok_auc],
                      paired = TRUE, alternative = "greater")  # H1: AUC(Boruta) > AUC(None)
    pt_acc <- t.test(w_acc[[bor_lab]][ok_acc], w_acc[[none_lab]][ok_acc],
                      paired = TRUE, alternative = "greater")  # H1: Acc(Boruta) > Acc(None)

    boruta_ttest_df <- data.frame(
      Metric           = c("AUC_Test", "Accuracy_Test"),
      N_pairs          = c(sum(ok_auc), sum(ok_acc)),
      Mean_Boruta      = c(mean(w_auc[[bor_lab]][ok_auc]),  mean(w_acc[[bor_lab]][ok_acc])),
      Mean_None        = c(mean(w_auc[[none_lab]][ok_auc]), mean(w_acc[[none_lab]][ok_acc])),
      Mean_Diff        = c(unname(pt_auc$estimate), unname(pt_acc$estimate)),
      t_statistic      = c(unname(pt_auc$statistic), unname(pt_acc$statistic)),
      df               = c(unname(pt_auc$parameter), unname(pt_acc$parameter)),
      p_value_one_tail = c(pt_auc$p.value, pt_acc$p.value),
      H1               = c("AUC higher with Boruta", "Accuracy higher with Boruta"),
      Significant_0.05 = c(pt_auc$p.value, pt_acc$p.value) < 0.05,
      stringsAsFactors = FALSE)

    # two-sided p: the one-tailed test only evaluates 'Boruta improves'; if Boruta is SIGNIFICANTLY
    # WORSE (one-tailed p close to 1), this two-sided p makes it visible.
    boruta_ttest_df$p_value_two_sided <- 2 * pt(-abs(boruta_ttest_df$t_statistic), boruta_ttest_df$df)
    write.xlsx(boruta_ttest_df, "Results/Excel/Boruta_Paired_Ttest.xlsx", overwrite = TRUE)
    cat("   OK Paired t-test exported to Results/Excel/Boruta_Paired_Ttest.xlsx\n")
    print(boruta_ttest_df)
  } else {
    boruta_ttest_df <- NULL
    cat("   Warning: The 2 Boruta levels ('boruta'/'none') could not be identified automatically; check results_df$Boruta manually.\n")
  }

  write.xlsx(list("Complete_Ranking"      = ranking, "Top_20" = head(ranking, 20),
                  "Summary_By_Model"      = ms, "Summary_By_Preprocessing" = ps,
                  "Summary_By_Boruta"     = bs),
             "Results/Excel/Ranking.xlsx", overwrite = TRUE)
  cat("OK Ranking exported\n")

  top10_plot <- ranking %>% head(10) %>%
    mutate(Combo = reorder(paste(Model, Preprocessing, Boruta, sep = " | "), AUC_Test)) %>%
    ggplot(aes(x = Combo, y = AUC_Test, fill = Model)) +
    geom_col() + geom_text(aes(label = round(AUC_Test, 4)), hjust = -0.1, size = 3.2) +
    scale_y_continuous(limits = c(0, 1.05)) + coord_flip() +
    labs(title = "Top 10 Best Combinations (AUC Test)", x = "", y = "AUC (Test set)") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  ggsave("Results/Other/Top10_Combinations.png", top10_plot, width = 13, height = 7, dpi = 300, bg = "white")
  cat("OK All plots generated\n")

} else { cat("\nError: No results generated\n") }

if (length(error_log) > 0) {
  write.xlsx(bind_rows(lapply(error_log, as.data.frame)),
             "Results/Excel/Error_Log.xlsx", overwrite = TRUE)
  cat("\nWarning: Errors:", length(error_log), "\n")
}

# ======================================
# 11. Final summary
# ======================================
cat("\n", rep("=", 70), "\n", sep = "")
cat("BINARY CLASSIFICATION PIPELINE COMPLETED\n")
cat(rep("=", 70), "\n", sep = "")

end_time <- Sys.time()
duration      <- difftime(end_time, start_time, units = "mins")

if (length(results_all) > 0) {
  cat("\n EXECUTION SUMMARY:\n")
  cat("  - Start:",    format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("  - End:",      format(end_time,    "%Y-%m-%d %H:%M:%S"), "\n")
  cat("  - Duration:", round(duration, 2), "minutes\n")
  cat("  - Samples after outlier removal:", nrow(X), "\n")
  cat("  - Outliers removed:", n_outliers, "(T2 AND Q)\n")
  cat("  - Classes:", paste(class_names, collapse = " vs "), "\n")
  cat("  - Positive class (AUC):", POS_CLASS, "\n")
  cat("  - Successful combinations:", nrow(results_df), "\n")
  best <- results_df %>% arrange(desc(AUC_Test), desc(Accuracy_Test), desc(Kappa_Test), N_vars) %>% head(1)
  cat("\n BEST COMBINATION (by AUC Test):\n")
  cat("  - Model:",         best$Model, "\n")
  cat("  - Preprocessing:", best$Preprocessing, "\n")
  cat("  - Boruta:",        best$Boruta, "\n")
  cat("  - AUC Test:",      round(best$AUC_Test, 4), "\n")
  cat("  - Accuracy Test:", round(best$Accuracy_Test, 4), "\n")
  cat("  - F1 Test:",       round(best$F1_Test, 4), "\n")
  cat("  - Kappa Test:",    round(best$Kappa_Test, 4), "\n")
}

cat("\n FILES GENERATED:\n")
cat("   Outliers/Raw_Spectra_Outliers.png\n")
cat("   Outliers/PCA_PC1_vs_PC2.png\n")
cat("   Outliers/PCA_PC2_vs_PC3.png\n")
cat("   Outliers/Influence_Plot_T2_vs_Q.png\n")
cat("   Heatmaps/Metrics_Model_Preprocessing.png\n")
cat("   Heatmaps/Metrics_Model_Boruta.png\n")
cat("   Heatmaps/Metrics_by_Model.png\n")
cat("   Heatmaps/Metrics_by_Preprocessing.png\n")
cat("   Heatmaps/AUC_Model_Preprocessing.png\n")
cat("   ConfusionMatrices/<Model>/  — Train & Test\n")
cat("   ROC/<Model>/               — ROC curve (Test)\n")
cat("   Other/Model_Comparison_AUC.png\n")
cat("   Other/Model_Comparison_Accuracy.png\n")
cat("   Other/Top10_Combinations.png\n")
cat("   Other/PCA_Metrics.png\n")
cat("   Excel/Complete_Summary.xlsx\n")
cat("   Excel/Ranking.xlsx\n")
cat("   Excel/PCA_Analysis.xlsx\n")

stopCluster(cl)
registerDoSEQ()

write.table(
  data.frame(Start    = format(start_time, "%Y-%m-%d %H:%M:%S"),
             End      = format(end_time,    "%Y-%m-%d %H:%M:%S"),
             Duration_min     = round(duration, 2),
             Samples_clean    = nrow(X),
             Outliers_removed = n_outliers,
             Criterion        = "T2 AND Q",
             Positive_class   = POS_CLASS,
             Combinations     = ifelse(length(results_all) > 0, nrow(results_df), 0),
             Errors           = length(error_log)),
  "Results/Execution_Info.txt", row.names = FALSE, quote = FALSE)

cat("\n All done! Check 'Results/'\n\n")

if (length(error_log) > 0) {
  print(bind_rows(lapply(error_log, as.data.frame)))
}




# ============================================================
# SECTION 12 — PDF REPORT  /  BINARY CLASSIFICATION PIPELINE
# Runs at the end of the pipeline, after cat("All done!")
# Requires: install.packages(c("grid","gridExtra","png"))
# ============================================================

library(grid)
library(gridExtra)
library(png)

# ── PALETTE ──────────────────────────────────────────────────
COL_DARK   <- "#1a1a2e"
COL_MID    <- "#16213e"
COL_LIGHT  <- "#0f3460"
COL_ACCENT <- "#4a90d9"
COL_SILVER <- "#a8c8e8"

gp_section <- gpar(fontsize = 13, fontface = "bold",  col = "white")
gp_sub     <- gpar(fontsize = 12, fontface = "bold",  col = COL_LIGHT)
gp_body    <- gpar(fontsize = 10, fontface = "plain", col = "#2c2c2c")
gp_caption <- gpar(fontsize =  8, fontface = "italic",col = "gray45")
gp_mono    <- gpar(fontsize =  9, fontface = "plain", col = "gray40",
                   fontfamily = "mono")

# ── HELPERS ──────────────────────────────────────────────────
insert_png <- function(path) {
  if (!file.exists(path))
    return(textGrob(paste0("[Missing: ", basename(path), "]"),
                    gp = gpar(col = "red", fontsize = 9)))
  rasterGrob(readPNG(path), interpolate = TRUE,
             width = unit(1,"npc"), height = unit(1,"npc"))
}

draw_page_number <- function(n) {
  grid.lines(x = c(0.05,0.95), y = c(0.028,0.028),
             gp = gpar(col = "gray80", lwd = 0.5))
  grid.text(paste("Page", n), x = unit(.5,"npc"), y = unit(.016,"npc"),
            gp = gpar(fontsize = 9, col = "gray50"))
}

draw_section_bar <- function(label, y_top = 0.955) {
  grid.rect(x = unit(.04,"npc"), y = unit(y_top,"npc"),
            width = unit(.92,"npc"), height = unit(.048,"npc"),
            just = c("left","top"), gp = gpar(fill = COL_MID, col = NA))
  grid.text(paste0("  ", label),
            x = unit(.05,"npc"), y = unit(y_top - .024,"npc"),
            just = c("left","center"), gp = gp_section)
}

draw_body <- function(txt, x = .06, y = .88, w = 108, lh = .034) {
  lines <- strwrap(txt, width = w)
  for (i in seq_along(lines))
    grid.text(lines[i], x = unit(x,"npc"),
              y = unit(y - (i-1)*lh,"npc"),
              just = c("left","top"), gp = gp_body)
  invisible(y - length(lines)*lh)
}

hrule <- function(y, col = "gray80")
  grid.lines(x = c(.05,.95), y = c(y,y), gp = gpar(col=col, lwd=.6))

# ── Diagram helpers ──────────────────────────────────────────
box_r <- function(x, y, w=.14, h=.055, fill=COL_ACCENT, label="", sz=8.5) {
  grid.roundrect(x=unit(x,"npc"), y=unit(y,"npc"),
                 width=unit(w,"npc"), height=unit(h,"npc"), r=unit(4,"pt"),
                 gp=gpar(fill=fill, col="white", lwd=1))
  grid.text(label, x=unit(x,"npc"), y=unit(y,"npc"),
            gp=gpar(fontsize=sz, fontface="bold", col="white"))
}
arrowh <- function(x0,x1,y,col=COL_ACCENT)
  grid.lines(x=c(x0,x1), y=c(y,y),
             arrow=arrow(length=unit(5,"pt"), type="closed"),
             gp=gpar(col=col, lwd=1.5, fill=col))

# ── Diagrams per algorithm ────────────────────────────────────
draw_diag_lda <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc", col="#dde6f0", lwd=.8))
  grid.text("LDA — Linear Discriminant Analysis",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  # two class clouds
  set.seed(1)
  c1x<-rnorm(18,.30,.07); c1y<-rnorm(18,.62,.08)
  c2x<-rnorm(18,.70,.07); c2y<-rnorm(18,.38,.08)
  grid.circle(x=c1x,y=c1y,r=unit(4,"pt"),
              gp=gpar(fill="#377EB8",col="white",lwd=.4,alpha=.85))
  grid.circle(x=c2x,y=c2y,r=unit(4,"pt"),
              gp=gpar(fill="#E41A1C",col="white",lwd=.4,alpha=.85))
  # decision boundary
  grid.lines(x=c(.22,.78),y=c(.78,.22),gp=gpar(col=COL_MID,lwd=2))
  grid.text("Decision\nboundary",x=.55,y=.70,
            gp=gpar(fontsize=7.5,col=COL_MID,fontface="bold"))
  grid.text("Class 1",x=.20,y=.82,gp=gpar(fontsize=8,col="#377EB8",fontface="bold"))
  grid.text("Class 2",x=.80,y=.20,gp=gpar(fontsize=8,col="#E41A1C",fontface="bold"))
  grid.text("Maximizes between-class /\nminimizes within-class variance",
            x=.50,y=.10,gp=gpar(fontsize=7,col="gray50"))
  popViewport()
}

draw_diag_plsda <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("PLS-DA — Latent Variable Classifier",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  box_r(.15,.68,w=.13,h=.10,fill=COL_LIGHT,label="X matrix\n(spectra)")
  box_r(.15,.44,w=.13,h=.08,fill="#c0392b",label="y dummy\n(0/1)")
  box_r(.45,.68,w=.16,h=.10,fill=COL_ACCENT,label="T scores\n(latent)")
  box_r(.45,.44,w=.16,h=.08,fill="#e67e22",label="U scores\n(class)")
  box_r(.75,.56,w=.14,h=.08,fill="#27ae60",label="Class\nprediction")
  arrowh(.22,.37,.68,COL_ACCENT); arrowh(.22,.37,.44,"#e67e22")
  arrowh(.53,.68,.68,COL_ACCENT); arrowh(.53,.68,.44,"#e67e22")
  grid.lines(x=c(.75,.75),y=c(.60,.68),gp=gpar(col=COL_ACCENT,lwd=1.2,lty="dashed"))
  grid.lines(x=c(.75,.75),y=c(.52,.44),gp=gpar(col="#e67e22",lwd=1.2,lty="dashed"))
  grid.text("ncomp tuned\nby 5-fold CV",x=.45,y=.22,
            gp=gpar(fontsize=7.5,col=COL_LIGHT,fontface="italic"))
  popViewport()
}

draw_diag_svm <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("SVM-RBF — Maximum Margin Classifier",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  set.seed(42)
  bx<-runif(14,.10,.48); by<-runif(14,.35,.78)
  rx<-runif(12,.52,.90); ry<-runif(12,.30,.75)
  grid.circle(x=bx,y=by,r=unit(4,"pt"),
              gp=gpar(fill=COL_ACCENT,col="white",lwd=.5,alpha=.8))
  grid.circle(x=rx,y=ry,r=unit(4,"pt"),
              gp=gpar(fill="#c0392b",col="white",lwd=.5,alpha=.8))
  grid.lines(x=c(.46,.54),y=c(.82,.18),gp=gpar(col=COL_MID,lwd=2))
  grid.lines(x=c(.38,.46),y=c(.82,.18),gp=gpar(col=COL_MID,lwd=1,lty="dashed"))
  grid.lines(x=c(.54,.62),y=c(.82,.18),gp=gpar(col=COL_MID,lwd=1,lty="dashed"))
  grid.text("C & \u03c3\nauto-tuned",x=.18,y=.16,
            gp=gpar(fontsize=7.5,col=COL_ACCENT))
  grid.text("RBF kernel:\nK(x,x')",x=.80,y=.16,
            gp=gpar(fontsize=7.5,col="#c0392b"))
  popViewport()
}

draw_diag_rf <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("Random Forest — Bootstrap Ensemble",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  tree_x<-c(.18,.50,.82)
  fills <-c("#2980b9","#27ae60","#8e44ad")
  for(i in 1:3){
    tx<-tree_x[i]
    grid.lines(x=c(tx,tx),y=c(.28,.42),gp=gpar(col=fills[i],lwd=2))
    grid.lines(x=c(tx-.10,tx,tx+.10),y=c(.58,.42,.58),gp=gpar(col=fills[i],lwd=1.5))
    grid.lines(x=c(tx-.06,tx-.10),y=c(.68,.58),gp=gpar(col=fills[i],lwd=1.2))
    grid.lines(x=c(tx+.06,tx+.10),y=c(.68,.58),gp=gpar(col=fills[i],lwd=1.2))
    grid.circle(x=tx,y=.72,r=unit(12,"pt"),gp=gpar(fill=fills[i],col=NA,alpha=.85))
    grid.text(paste0("Tree ",i),x=tx,y=.20,
              gp=gpar(fontsize=7.5,col=fills[i],fontface="bold"))
    grid.text("Bootstrap",x=tx,y=.12,gp=gpar(fontsize=6.5,col="gray50"))
  }
  grid.lines(x=c(.18,.50,.82),y=c(.28,.22,.28),
             gp=gpar(col=COL_DARK,lwd=1.2,lty="dashed"))
  box_r(.50,.10,w=.20,h=.07,fill=COL_DARK,label="Majority vote",sz=7.5)
  popViewport()
}

draw_diag_dt <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("Decision Tree (CART) — Recursive Partitioning",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  # root
  box_r(.50,.80,w=.22,h=.07,fill=COL_MID,label="Root node\n(best split)",sz=8)
  # level 1
  box_r(.28,.60,w=.20,h=.07,fill=COL_ACCENT,label="Node A",sz=8)
  box_r(.72,.60,w=.20,h=.07,fill=COL_ACCENT,label="Node B",sz=8)
  # level 2
  box_r(.16,.38,w=.18,h=.07,fill="#27ae60",label="Leaf\nClass 1",sz=7.5)
  box_r(.38,.38,w=.18,h=.07,fill="#c0392b",label="Leaf\nClass 2",sz=7.5)
  box_r(.62,.38,w=.18,h=.07,fill="#27ae60",label="Leaf\nClass 1",sz=7.5)
  box_r(.84,.38,w=.18,h=.07,fill="#c0392b",label="Leaf\nClass 2",sz=7.5)
  # lines
  grid.lines(x=c(.50,.28),y=c(.765,.635),gp=gpar(col=COL_MID,lwd=1))
  grid.lines(x=c(.50,.72),y=c(.765,.635),gp=gpar(col=COL_MID,lwd=1))
  grid.lines(x=c(.28,.16),y=c(.565,.415),gp=gpar(col=COL_ACCENT,lwd=1))
  grid.lines(x=c(.28,.38),y=c(.565,.415),gp=gpar(col=COL_ACCENT,lwd=1))
  grid.lines(x=c(.72,.62),y=c(.565,.415),gp=gpar(col=COL_ACCENT,lwd=1))
  grid.lines(x=c(.72,.84),y=c(.565,.415),gp=gpar(col=COL_ACCENT,lwd=1))
  grid.text("cp (complexity) tuned by 5-fold CV",
            x=.50,y=.18,gp=gpar(fontsize=7.5,col="gray50",fontface="italic"))
  popViewport()
}

draw_diag_knn <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("k-Nearest Neighbours (kNN)",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  set.seed(9)
  c1x<-c(.20,.25,.18,.30,.22); c1y<-c(.65,.72,.60,.68,.75)
  c2x<-c(.70,.75,.68,.80,.72); c2y<-c(.40,.35,.50,.42,.38)
  # query point
  qx<-.45; qy<-.57
  grid.circle(x=c1x,y=c1y,r=unit(5,"pt"),gp=gpar(fill="#377EB8",col="white",lwd=.5))
  grid.circle(x=c2x,y=c2y,r=unit(5,"pt"),gp=gpar(fill="#E41A1C",col="white",lwd=.5))
  # k=3 neighbours
  nx<-c(.38,.42,.50); ny<-c(.65,.48,.62)
  nc<-c("#377EB8","#377EB8","#E41A1C")
  grid.circle(x=nx,y=ny,r=unit(5,"pt"),gp=gpar(fill=nc,col="white",lwd=.5))
  # query
  grid.circle(x=qx,y=qy,r=unit(7,"pt"),gp=gpar(fill="gold",col=COL_DARK,lwd=1.5))
  grid.text("?",x=qx,y=qy,gp=gpar(fontsize=9,fontface="bold",col=COL_DARK))
  # distances
  for(i in 1:3)
    grid.lines(x=c(qx,nx[i]),y=c(qy,ny[i]),
               gp=gpar(col="gray50",lwd=.8,lty="dashed"))
  # circle
  grid.circle(x=qx,y=qy,r=unit(32,"pt"),
              gp=gpar(col=COL_MID,fill=NA,lwd=1,lty="dashed"))
  grid.text("k = 3, 5, 7, 9, 11, 13, 15\n(majority vote within k neighbours)",
            x=.50,y=.12,gp=gpar(fontsize=7.5,col="gray50"))
  popViewport()
}

draw_diag_nb <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("Naive Bayes — P(Class | X) via Bayes Theorem",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  # Bayes formula illustration
  grid.text("P(C | X) \u221d P(X | C) \u00d7 P(C)",
            x=.50,y=.75,gp=gpar(fontsize=12,fontface="bold",col=COL_MID))
  hrule_y<-.67
  grid.lines(x=c(.10,.90),y=c(hrule_y,hrule_y),gp=gpar(col="gray70",lwd=.6))
  # three boxes
  box_r(.22,.52,w=.24,h=.10,fill=COL_ACCENT,label="Likelihood\nP(X | C)",sz=8)
  box_r(.50,.52,w=.22,h=.10,fill="#e67e22",label="Prior\nP(C)",sz=8)
  box_r(.78,.52,w=.24,h=.10,fill="#27ae60",label="Posterior\nP(C | X)",sz=8)
  grid.text("\u00d7",x=.37,y=.52,gp=gpar(fontsize=14,col=COL_DARK,fontface="bold"))
  grid.text("\u2192",x=.64,y=.52,gp=gpar(fontsize=14,col=COL_DARK,fontface="bold"))
  grid.text("Kernel density used (usekernel=TRUE)\nLaplace smoothing: 0, 1",
            x=.50,y=.24,gp=gpar(fontsize=7.5,col="gray50"))
  grid.text("\"Naive\" = features assumed independent given class",
            x=.50,y=.12,gp=gpar(fontsize=7,col="gray50",fontface="italic"))
  popViewport()
}

draw_diag_xgb <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("XGBoost — Sequential Residual Boosting",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  fills<-c("#e74c3c","#e67e22","#f1c40f","#2ecc71")
  xs<-seq(.12,.82,length.out=4)
  for(i in 1:4){
    h<-.06+.03*(4-i)
    grid.rect(x=unit(xs[i],"npc"),y=unit(.65,"npc"),
              width=unit(.10,"npc"),height=unit(h*3,"npc"),
              gp=gpar(fill=fills[i],col=NA,alpha=.85))
    grid.text(paste0("T",i),x=xs[i],y=.65,
              gp=gpar(fontsize=8,fontface="bold",col="white"))
    grid.text(paste0("r",i),x=xs[i],y=.44,
              gp=gpar(fontsize=6.5,col=fills[i]))
    if(i<4) arrowh(xs[i]+.05,xs[i+1]-.05,.65,COL_DARK)
  }
  grid.lines(x=c(.12,.82),y=c(.30,.30),gp=gpar(col=COL_DARK,lwd=1,lty="dashed"))
  box_r(.50,.18,w=.30,h=.07,fill=COL_DARK,
        label="logistic(T1+T2+T3+T4)",sz=7.5)
  grid.text("nthread=1 | direct xgb.train (no caret wrapper)",
            x=.50,y=.06,gp=gpar(fontsize=7,col="gray50",fontface="italic"))
  popViewport()
}

# ══════════════════════════════════════════════════════════════
# BORUTA LOG — ensure it exists even if empty
# ══════════════════════════════════════════════════════════════
# The pipeline should generate boruta_vars_log during the main loop.
# If for any reason it does not exist, it is initialized empty.
if (!exists("boruta_vars_log")) boruta_vars_log <- list()

# ══════════════════════════════════════════════════════════════
# GENERATE SPECTRA + BORUTA PNGs (before opening the PDF)
# One PNG per preprocessing: raw spectrum (X_train, post-outlier,
# no preprocessing) with translucent red lines at the
# wavenumbers confirmed by Boruta.
# ══════════════════════════════════════════════════════════════
cat(">>> Generating Boruta spectra overlay plots (raw spectra)...\n")
dir.create("Results/Spectra/Boruta_Overlays", showWarnings=FALSE, recursive=TRUE)

if (!exists("boruta_vars_log")) boruta_vars_log <- list()

# Function to build the spectra plot with Boruta marks
make_spectra_boruta_plot_clf <- function(X_matrix, wl_vec,
                                         sel_vars  = NULL,
                                         title_str = "",
                                         y_classes = NULL) {
  Xm  <- as.matrix(X_matrix)
  nc  <- ncol(Xm); nr <- nrow(Xm)
  wl  <- if(length(wl_vec)==nc) wl_vec else seq_len(nc)
  
  df <- data.frame(
    Wavenumber = rep(wl, each=nr),
    Absorbance = as.vector(Xm),
    Sample     = rep(seq_len(nr), times=nc)
  )
  df <- df[complete.cases(df),]
  
  # Color by class if available
  use_class <- !is.null(y_classes) && length(y_classes)==nr
  if(use_class) df$Class <- rep(as.character(y_classes), times=nc)
  
  cls_pal <- c("#377EB8","#E41A1C","#4DAF4A","#984EA3","#FF7F00")
  
  p <- if(use_class){
    ggplot(df, aes(x=Wavenumber, y=Absorbance, group=Sample, color=Class)) +
      scale_color_manual(values=setNames(cls_pal[seq_along(unique(df$Class))],
                                         unique(df$Class)),
                         name="Class")
  } else {
    ggplot(df, aes(x=Wavenumber, y=Absorbance, group=Sample))
  }
  
  p <- p +
    geom_line(alpha=.28, linewidth=.35) +
    labs(title    = title_str,
         x        = expression("Wavenumber (cm"^{-1}*")"),
         y        = "Absorbance") +
    theme_minimal(base_size=12) +
    theme(plot.title       = element_text(face="bold", hjust=.5, size=13),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color="gray70", fill=NA),
          legend.position  = if(use_class) "right" else "none")
  
  # Translucent red vertical lines
  if(!is.null(sel_vars) && length(sel_vars)>0){
    sw  <- suppressWarnings(parse_wn(sel_vars))
    sw  <- sw[!is.na(sw) & sw %in% wl]
    if(length(sw)>0){
      p <- p +
        geom_vline(data=data.frame(xintercept=sw),
                   aes(xintercept=xintercept),
                   color="red", alpha=.35, linewidth=.9,
                   inherit.aes=FALSE) +
        annotate("text",
                 x=min(wl,na.rm=TRUE),
                 y=max(df$Absorbance,na.rm=TRUE),
                 label=paste0(length(sw)," variables selected (Boruta)"),
                 hjust=0, vjust=1, color="red", size=3.8, fontface="bold")
    }
  }
  p
}

# Training set classes used for coloring
y_cls_train <- if(exists("y_train")) as.character(y_train) else NULL

# Overview PNG without marks
raw_spectra_clf_path <- "Results/Spectra/Boruta_Overlays/_raw_all_samples.png"
if(exists("X_train") && exists("wavelengths")){
  p_ov <- make_spectra_boruta_plot_clf(
    X_matrix  = X_train,
    wl_vec    = wavelengths,
    sel_vars  = NULL,
    title_str = paste0("Raw Spectra — Training set (",nrow(X_train),
                       " samples). Colour = class."),
    y_classes = y_cls_train
  )
  tryCatch(ggsave(raw_spectra_clf_path, p_ov, width=12, height=5, dpi=150,bg="white"),
           error=function(e) NULL)
  cat("   OK Raw overview saved\n")
}

# One PNG per preprocessing: raw + Boruta
boruta_png_map_clf <- list()

for(sc in scatter_opts){
  for(dv in deriv_opts){
    pp_nm_loc <- make_pp_name(sc, dv)
    
    bor_key <- paste(pp_nm_loc, "boruta", sep="__")
    
    sel_vars_loc <- NULL
    if(bor_key %in% names(boruta_vars_log)){
      bvl <- boruta_vars_log[[bor_key]]
      if(nchar(bvl$Variables)>0 && !isTRUE(bvl$use_all))
        sel_vars_loc <- trimws(strsplit(bvl$Variables,",")[[1]])
    }
    
    path_rb <- paste0("Results/Spectra/Boruta_Overlays/",
                      gsub(" ","_",pp_nm_loc),"_raw_boruta.png")
    
    p_rb <- make_spectra_boruta_plot_clf(
      X_matrix  = X_train,
      wl_vec    = wavelengths,
      sel_vars  = sel_vars_loc,
      title_str = paste0("Raw spectra + Boruta selection — ", pp_nm_loc),
      y_classes = y_cls_train
    )
    tryCatch(ggsave(path_rb, p_rb, width=12, height=5, dpi=150,bg="white"),
             error=function(e) NULL)
    
    boruta_png_map_clf[[pp_nm_loc]] <- list(
      raw_boruta = path_rb,
      n_vars     = if(!is.null(sel_vars_loc)) length(sel_vars_loc) else NA
    )
    cat("   OK", pp_nm_loc,
        if(!is.null(sel_vars_loc)) paste0("(",length(sel_vars_loc)," vars)") else "(all)",
        "\n")
  }
}
cat(">>> Done.\n\n")

# ══════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════
pdf_path <- "Results/Scan_Classification_Report.pdf"
pdf(pdf_path, width=11, height=8.5, paper="USr")
pg <- 0L

# ══════════════════════════════════════════════════════════════
# P1 — COVER PAGE
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L

grid.rect(x=0, y=1, width=1, height=.52, just=c("left","top"),
          gp=gpar(fill=COL_DARK, col=NA))

grid.text("FT-MIR Spectroscopy",
          x=.5, y=.90,
          gp=gpar(fontsize=26,fontface="bold",col="white"))
grid.text("Pipeline for Binary Classification Screening",
          x=.5, y=.81,
          gp=gpar(fontsize=19,fontface="plain",col=COL_SILVER))
grid.text("Machine Learning-Based Qualitative Analysis Report",
          x=.5, y=.73,
          gp=gpar(fontsize=12,fontface="italic",col="#c8d8e8"))
grid.lines(x=c(.10,.90),y=c(.685,.685),gp=gpar(col=COL_SILVER,lwd=1.2))

# ── Metadata block ─────────────────────────────────────────
grid.rect(x=unit(.06,"npc"), y=unit(.08,"npc"),
          width=unit(.52,"npc"), height=unit(.37,"npc"),
          just=c("left","bottom"),
          gp=gpar(fill="white",col=COL_MID,lwd=1.2))
grid.rect(x=unit(.06,"npc"), y=unit(.45,"npc"),
          width=unit(.52,"npc"), height=unit(.028,"npc"),
          just=c("left","top"), gp=gpar(fill=COL_MID,col=NA))
grid.text("Study Information",
          x=unit(.32,"npc"), y=unit(.436,"npc"),
          gp=gpar(fontsize=10,fontface="bold",col="white"))

REPORT_MATRIX     <- ""
REPORT_ANALYTE    <- sheet_name   # filled in automatically from the analyzed sheet
REPORT_INSTRUMENT <- ""
REPORT_N_SAMPLES  <- if(exists("X_raw")) as.character(nrow(X_raw)) else ""
REPORT_CLASSES    <- if(exists("class_names")) paste(class_names, collapse=" vs ") else ""

fields <- list(
  list(label="Matrix:",                 val=REPORT_MATRIX),
  list(label="Analyte / Target:",       val=REPORT_ANALYTE),
  list(label="Instrument:",             val=REPORT_INSTRUMENT),
  list(label="Number of total samples:",val=REPORT_N_SAMPLES),
  list(label="Classes:",                val=REPORT_CLASSES)
)
for(i in seq_along(fields)){
  yy <- .395 - (i-1)*.065
  grid.text(fields[[i]]$label,
            x=unit(.10,"npc"), y=unit(yy,"npc"),
            just=c("left","center"),
            gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
  grid.lines(x=c(.30,.555),y=c(yy-.012,yy-.012),
             gp=gpar(col="gray70",lwd=.7,lty="dotted"))
  if(nchar(fields[[i]]$val)>0)
    grid.text(fields[[i]]$val,
              x=unit(.31,"npc"), y=unit(yy,"npc"),
              just=c("left","center"),
              gp=gpar(fontsize=10,col="#333333"))
}

grid.text(paste("Generated:",format(Sys.time(),"%Y-%m-%d %H:%M")),
          x=unit(.94,"npc"), y=unit(.12,"npc"),
          just=c("right","center"),gp=gpar(fontsize=9,col="gray55"))
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P2 — TABLE OF CONTENTS
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
grid.rect(x=0,y=1,width=1,height=.12,just=c("left","top"),
          gp=gpar(fill=COL_MID,col=NA))
grid.text("Table of Contents",x=.5,y=.94,
          gp=gpar(fontsize=18,fontface="bold",col="white"))

toc <- list(
  list(n="1",  t="Outlier Detection & Removal",                          p=3),
  list(n="2",  t="Spectral Preprocessing Methods",                       p=5),
  list(n="3",  t="Data Split & Cross-Validation Strategy",               p=6),
  list(n="4",  t="Machine Learning Algorithms & Hyperparameters",        p=7),
  list(n="4",  t="  \u2514 Hyperparameter Summary Table (all models)",   p=15),
  list(n="5",  t="Boruta Variable Selection — Summary & Variable Lists", p=16),
  list(n="6",  t="Boruta Variables Highlighted on Raw Spectra",          p=18),
  list(n="7",  t="Performance Results — Summary Tables",                 p=29),
  list(n="8",  t="Heatmaps — Model vs Preprocessing",                    p=31),
  list(n="9",  t="Confusion Matrices — by Model",                        p=39),
  list(n="10", t="ROC Curves — by Model",                                p=47)
)
yy <- .83
for(e in toc){
  grid.text(paste0(e$n,".   ",e$t),
            x=unit(.08,"npc"),y=unit(yy,"npc"),
            just=c("left","center"),gp=gpar(fontsize=10.5,col=COL_DARK))
  grid.text(as.character(e$p),
            x=unit(.92,"npc"),y=unit(yy,"npc"),
            just=c("right","center"),gp=gpar(fontsize=10.5,col="gray50"))
  grid.lines(x=c(.42,.88),y=c(yy-.008,yy-.008),
             gp=gpar(col="gray82",lwd=.5,lty="dotted"))
  yy <- yy - .06
}
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P3 — OUTLIERS: text
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("1. Outlier Detection & Removal")

t1 <- paste0(
  "Outlier detection was performed on the raw spectral matrix (",
  nrow(X_raw)," samples \u00d7 ",ncol(X_raw)," variables) using a two-criterion ",
  "PCA-based approach. Data were mean-centered; the number of PCs retained was the ",
  "minimum needed to explain \u226595% of cumulative spectral variance (",
  n_pcs," PCs, ",round(var_cum[n_pcs]*100,1),"% variance).")
y_cur <- draw_body(t1, y=.87)

t2 <- paste0(
  "Two criteria applied simultaneously: (1) Hotelling T\u00b2 compared against a ",
  "chi-squared threshold at 99% confidence (T\u00b2 limit = ",round(T2_lim,2),
  "). (2) Q residuals exceeding mean + 3 SD (Q limit = ",round(Q_lim,4),
  "). Only samples exceeding BOTH thresholds simultaneously (\"severe\" outliers) were removed. ",
  n_outliers," sample(s) removed \u2192 ",sum(!flag_outlier)," retained.")
y_cur <- draw_body(t2, y=y_cur-.025)

# summary box
grid.rect(x=unit(.05,"npc"),y=unit(y_cur-.075,"npc"),
          width=unit(.90,"npc"),height=unit(.048,"npc"),
          just=c("left","bottom"),
          gp=gpar(fill="#eef4fb",col="#b0c8e0",lwd=.8))
grid.text(
  paste0("T\u00b2 flagged: ",sum(flag_T2),
         "   |   Q flagged: ",sum(flag_Q),
         "   |   Removed (T\u00b2 AND Q): ",n_outliers," of ",nrow(X_raw),
         "   |   Classes: ",paste(class_names,collapse=" vs ")),
  x=.5, y=y_cur-.052,
  gp=gpar(fontsize=9.5,fontface="bold",col=COL_MID))

img_raw <- insert_png("Results/Outliers/Raw_Spectra_Outliers.png")
grid.draw(editGrob(img_raw, vp=viewport(x=.5,y=.285,width=.90,height=.41)))
grid.text("Figure 1. Raw FTMIR spectra. Each class coloured separately. Outliers = dashed black.",
          x=.5,y=.065,gp=gp_caption)
draw_page_number(pg)

# ── P4 PCA biplots ────────────────────────────────────────────
grid.newpage(); pg <- pg + 1L
draw_section_bar("1. Outlier Detection & Removal (cont.)")

img_pc12 <- insert_png("Results/Outliers/PCA_PC1_vs_PC2.png")
img_pc23 <- insert_png("Results/Outliers/PCA_PC2_vs_PC3.png")
img_infl <- insert_png("Results/Outliers/Influence_Plot_T2_vs_Q.png")

grid.draw(editGrob(img_pc12,vp=viewport(x=.25,y=.67,width=.46,height=.44)))
grid.draw(editGrob(img_pc23,vp=viewport(x=.75,y=.67,width=.46,height=.44)))
grid.text("Figure 2. PCA PC1 vs PC2. Classes coloured.",x=.25,y=.44,gp=gp_caption)
grid.text("Figure 3. PCA PC2 vs PC3.",x=.75,y=.44,gp=gp_caption)
hrule(.425)
grid.draw(editGrob(img_infl,vp=viewport(x=.5,y=.245,width=.55,height=.36)))
grid.text("Figure 4. Influence plot: T\u00b2 vs Q residuals. Dashed lines = decision thresholds.",
          x=.5,y=.055,gp=gp_caption)
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P5 — PREPROCESSING
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("2. Spectral Preprocessing Methods")

preproc_items <- list(
  list(name="Mean Centering  (applied to all combinations)",
       desc=paste0("Each spectral variable centered by subtracting the training-set column mean. ",
                   "Same mean applied to test set to prevent data leakage. Prerequisite for all linear methods.")),
  list(name="Savitzky-Golay smoothing / derivatives  (0th, 1st, 2nd | poly=3, window=11)",
       desc=paste0("Computed on the complete, continuous spectrum (545 variables) before the Patz windows are cropped. ",
                   "0th = smoothed spectrum", if (SMOOTH_ALWAYS) " (applied to every combination, including those without derivative). " else " (NOT applied in this run: combinations without derivative are raw spectra). ",
                   "1st derivative removes additive baseline offsets. ",
                   "2nd derivative removes constant and linear baselines, enhances spectral resolution. ",
                   if (SCATTER_FIRST) "Order of operations: scatter correction (SNV/MSC) -> Savitzky-Golay -> Patz windows -> mean centering."
                   else "Order of operations: Savitzky-Golay -> scatter correction (SNV/MSC) -> Patz windows -> mean centering.")),
  list(name="Standard Normal Variate (SNV)",
       desc=paste0("Each spectrum scaled to zero mean and unit variance independently. ",
                   "Corrects multiplicative scatter and path-length differences. No reference spectrum required.")),
  list(name="Multiplicative Scatter Correction (MSC)",
       desc=paste0("Linear regression of each spectrum against the training mean spectrum; additive and ",
                   "multiplicative scatter estimated and removed. Reference computed once from training set only."))
)

yy <- .865
for(pp_item in preproc_items){
  grid.text(paste0("\u25B6  ",pp_item$name),
            x=unit(.06,"npc"),y=unit(yy,"npc"),
            just=c("left","center"),gp=gp_sub)
  yy <- draw_body(pp_item$desc, y=yy-.032) - .022
}
hrule(yy-.01)
grid.text(paste0("Total combinations: ",
                 length(scatter_opts)*length(deriv_opts),
                 "   (Scatter: none / SNV / MSC  \u00d7  Derivative: 0 / 1st / 2nd)"),
          x=.5, y=yy-.040,
          gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P6 — DATA SPLIT & CV
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("3. Data Split & Cross-Validation Strategy")

n_tr <- if(exists("X_train")) nrow(X_train) else "N/A"
n_te <- if(exists("X_test"))  nrow(X_test)  else "N/A"

y_cur <- draw_body(paste0(
  "After outlier removal, the ",sum(!flag_outlier),"-sample dataset was partitioned ",
  "into training and test sets using a stratified 70/30 split (createDataPartition, p=0.70, ",
  "set.seed=1234). Stratification on the class label ensures balanced class proportions in both ",
  "subsets, which is critical when class imbalance exists."), y=.87)

y_cur <- draw_body(paste0(
  "All caret models used 5-fold cross-validation (trainControl, method='cv', number=5) with ",
  "classProbs=TRUE and twoClassSummary for ROC-based tuning. XGBoost was trained directly via ",
  "xgb.train() with manual 5-fold CV and nthread=1 to avoid the parallelism conflict. ",
  "The test set was held out entirely and used only for final evaluation."),
  y=y_cur-.025)

split_df <- data.frame(
  Parameter=c("Samples after outlier removal","Training set","Test set",
              "Split ratio","Stratification","Random seed",
              "Cross-validation","Tuning metric (caret)",
              "Positive class (ROC/AUC)",
              "Train class distribution","Test class distribution"),
  Value=c(sum(!flag_outlier), n_tr, n_te, "70 / 30", "Yes — by class label", "1234",
          "5-fold CV (training set only)", "ROC (AUC)",
          if(exists("POS_CLASS")) POS_CLASS else "N/A",
          if(exists("y_train")) paste(names(table(y_train)),table(y_train),sep="=",collapse=" / ") else "N/A",
          if(exists("y_test"))  paste(names(table(y_test)), table(y_test), sep="=",collapse=" / ") else "N/A"),
  stringsAsFactors=FALSE)

tbl_s <- tableGrob(split_df, rows=NULL,
                   theme=ttheme_minimal(base_size=9.5,
                                        core   =list(fg_params=list(col=COL_DARK,hjust=0,x=.03)),
                                        colhead=list(fg_params=list(fontface="bold",col=COL_MID,hjust=0,x=.03))))
grid.draw(editGrob(tbl_s,
                   vp=viewport(x=.45,y=max(y_cur-.22,.25),width=.84,height=.40)))
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# SECTION 4 — ALGORITHMS (one page each)
# ══════════════════════════════════════════════════════════════
algo_list <- list(
  list(
    name="Linear Discriminant Analysis (LDA)",
    caret="method = 'lda',  preProcess = c('center','scale')",
    desc=paste0(
      "LDA finds the linear combination of features that maximises the ratio of between-class ",
      "to within-class variance (Fisher criterion). It projects the high-dimensional spectral ",
      "space onto a subspace where classes are maximally separated, then assigns each sample to ",
      "the nearest class centroid in that space. LDA assumes multivariate normality and equal ",
      "covariance matrices across classes. Despite these assumptions, it often performs well on ",
      "spectral data due to the high signal-to-noise ratio of FT-MIR spectra."),
    hyper=data.frame(
      Hyperparameter=c("No free hyperparameters"),
      Range=c("All parameters estimated analytically from training data"),
      stringsAsFactors=FALSE),
    diag=draw_diag_lda
  ),
  list(
    name="PLS-DA — Partial Least Squares Discriminant Analysis",
    caret="method = 'pls',  preProcess = c('center','scale')",
    desc=paste0(
      "PLS-DA encodes the class label as a binary dummy variable (0/1) and applies PLS regression. ",
      "Latent variables (components) are extracted to maximise covariance between the spectral ",
      "matrix X and the class indicator y. Samples are classified by thresholding the continuous ",
      "PLS prediction at 0.5. This makes PLS-DA the standard chemometric classification method, ",
      "especially effective when p >> n as in spectral datasets."),
    hyper=data.frame(
      Hyperparameter=c("ncomp (latent variables)"),
      Range=c("1 to 15  (tuneLength = 15)"),
      stringsAsFactors=FALSE),
    diag=draw_diag_plsda
  ),
  list(
    name="Support Vector Machine — Radial Basis (SVM-RBF)",
    caret="method = 'svmRadial',  preProcess = c('center','scale')",
    desc=paste0(
      "SVM-RBF maps inputs to a high-dimensional feature space via K(x,x') = exp(\u2212\u03c3||x\u2212x'||\u00b2) ",
      "and finds a maximum-margin hyperplane separating the classes. The cost C controls the ",
      "penalty for misclassifications; larger C = smaller margin = risk of overfitting. ",
      "\u03c3 (RBF width) is auto-estimated via sigest(). SVM produces class probabilities via ",
      "Platt scaling (classProbs=TRUE), enabling AUC computation."),
    hyper=data.frame(
      Hyperparameter=c("C (cost)","sigma (\u03c3, RBF width)"),
      Range=c("Auto-estimated (tuneLength = 8)","Auto via sigest() heuristic"),
      stringsAsFactors=FALSE),
    diag=draw_diag_svm
  ),
  list(
    name="Random Forest (RF)",
    caret="method = 'rf',  ntree = 300",
    desc=paste0(
      "Random Forest builds 300 CART trees, each on an independent bootstrap sample. ",
      "At each node only mtry randomly chosen predictors are considered, decorrelating the trees. ",
      "Class probabilities are estimated as the fraction of trees voting for each class, ",
      "enabling AUC computation. RF also provides per-variable mean-decrease-in-accuracy ",
      "importance scores, useful for identifying informative spectral bands."),
    hyper=data.frame(
      Hyperparameter=c("mtry (predictors per split)","ntree"),
      Range=c("Auto: 6 candidate values  (tuneLength = 6)","Fixed = 300"),
      stringsAsFactors=FALSE),
    diag=draw_diag_rf
  ),
  list(
    name="Decision Tree (CART)",
    caret="method = 'rpart'",
    desc=paste0(
      "A single CART tree recursively partitions the feature space by selecting the split that ",
      "maximises the Gini impurity reduction at each node. The complexity parameter cp ",
      "(cost-complexity pruning) controls tree depth: larger cp = simpler tree. ",
      "CART is highly interpretable but prone to overfitting; it is included as a baseline ",
      "and for comparison with the ensemble methods RF and XGBoost."),
    hyper=data.frame(
      Hyperparameter=c("cp (complexity parameter)"),
      Range=c("10^(\u22124) to 10^(\u22121),  20 log-spaced values"),
      stringsAsFactors=FALSE),
    diag=draw_diag_dt
  ),
  list(
    name="k-Nearest Neighbours (kNN)",
    caret="method = 'knn',  preProcess = c('center','scale')",
    desc=paste0(
      "kNN classifies each sample by majority vote among its k nearest training neighbours in ",
      "the (preprocessed) feature space, using Euclidean distance. No explicit model is fitted; ",
      "the entire training set is retained as the 'model'. kNN is sensitive to irrelevant or ",
      "noisy features (curse of dimensionality), which is why Boruta variable selection is ",
      "particularly valuable for this method with high-dimensional spectral data."),
    hyper=data.frame(
      Hyperparameter=c("k (number of neighbours)"),
      Range=c("3, 5, 7, 9, 11, 13, 15  [7 values]"),
      stringsAsFactors=FALSE),
    diag=draw_diag_knn
  ),
  list(
    name="Naive Bayes (NB)",
    caret="method = 'naive_bayes'",
    desc=paste0(
      "Naive Bayes applies Bayes' theorem assuming conditional independence of features given ",
      "the class. Kernel density estimation (usekernel=TRUE) is used instead of Gaussian ",
      "assumptions, making it more flexible for non-normal spectral data. Laplace smoothing ",
      "prevents zero-probability issues. Despite the 'naive' independence assumption, NB often ",
      "performs surprisingly well on high-dimensional spectral data."),
    hyper=data.frame(
      Hyperparameter=c("laplace","usekernel","adjust (bandwidth)"),
      Range=c("0, 1","TRUE (kernel density estimation)","0.5, 1.0, 1.5"),
      stringsAsFactors=FALSE),
    diag=draw_diag_nb
  ),
  list(
    name="XGBoost — eXtreme Gradient Boosting (direct xgb.train)",
    caret="xgb.train() direct  |  objective = 'binary:logistic'",
    desc=paste0(
      "XGBoost builds trees sequentially, each fitting the negative gradient (pseudo-residuals) ",
      "of the binary cross-entropy loss. Predictions are passed through a sigmoid to produce ",
      "class probabilities. L1/L2 regularization on leaf weights, min_child_weight=1, and ",
      "feature/row subsampling prevent overfitting. Implementation uses xgb.train() directly ",
      "(not via caret) with nthread=1 to avoid the OpenMP/fork conflict that causes ",
      "'Error: Stopping' when caret's xgbTree runs inside doParallel workers."),
    hyper=data.frame(
      Hyperparameter=c("nrounds","max_depth","eta","gamma",
                       "colsample_bytree","min_child_weight","subsample"),
      Range=c("50, 100, 150","2, 3, 4","0.05, 0.10","0, 0.1",
              "Fixed = 0.80","Fixed = 1","Fixed = 0.80"),
      stringsAsFactors=FALSE),
    diag=draw_diag_xgb
  )
)

for(al in algo_list){
  grid.newpage(); pg <- pg + 1L
  draw_section_bar(paste0("4. Algorithm: ",al$name))
  grid.text(paste0("caret / API: ",al$caret),
            x=unit(.06,"npc"),y=unit(.875,"npc"),
            just=c("left","center"),gp=gp_mono)
  hrule(.862)
  
  lines_d <- strwrap(al$desc, width=62)
  for(i in seq_along(lines_d))
    grid.text(lines_d[i],x=unit(.06,"npc"),
              y=unit(.835-(i-1)*.032,"npc"),
              just=c("left","top"),gp=gp_body)
  text_bot <- .835 - length(lines_d)*.032
  
  # diagram right
  al$diag(viewport(x=.77,y=.67,width=.43,height=.38))
  
  hyp_y <- min(text_bot - .035, .44)
  grid.text("Hyperparameters evaluated:",
            x=unit(.06,"npc"),y=unit(hyp_y,"npc"),
            just=c("left","center"),gp=gp_sub)
  
  tbl_h <- tableGrob(al$hyper, rows=NULL,
                     theme=ttheme_minimal(base_size=9.5,
                                          core   =list(fg_params=list(col=COL_DARK,hjust=0,x=.02)),
                                          colhead=list(fg_params=list(fontface="bold",col=COL_MID,hjust=0,x=.02))))
  rh <- .052; th <- nrow(al$hyper)*rh + .065
  grid.draw(editGrob(tbl_h,
                     vp=viewport(x=.50,y=hyp_y-th/2-.03,width=.88,height=th)))
  draw_page_number(pg)
}

# ── Hyperparameter summary table — all models ─────────────────
grid.newpage(); pg <- pg + 1L
draw_section_bar("4. Hyperparameter Summary — All Models")

draw_body(paste0(
  "The table below consolidates all hyperparameters evaluated across every classification model. ",
  "Fixed values were held constant; tuned values were explored via 5-fold CV on the training set ",
  "(metric: ROC/AUC). XGBoost used a full grid search with manual 5-fold CV via xgb.train()."),
  y=.875)

hyper_summary_clf <- data.frame(
  Model = c(
    "LDA",
    "PLS-DA","PLS-DA",
    "SVM-RBF","SVM-RBF",
    "RF","RF",
    "CART",
    "kNN",
    "Naive Bayes","Naive Bayes","Naive Bayes",
    "XGBoost","XGBoost","XGBoost","XGBoost","XGBoost","XGBoost","XGBoost"
  ),
  Hyperparameter = c(
    "(analytical — no tuning)",
    "ncomp","[method]",
    "C (cost)","sigma (\u03c3)",
    "mtry","ntree",
    "cp (complexity)",
    "k (neighbours)",
    "laplace","usekernel","adjust (bandwidth)",
    "nrounds","max_depth","eta","gamma","colsample_bytree","min_child_weight","subsample"
  ),
  Values_Range = c(
    "All parameters estimated from training data",
    "1, 2, 3 \u2026 15  (15 values)","pls",
    "8 values auto-estimated from data","auto via sigest() heuristic",
    "6 values auto-selected (tuneLength=6)","Fixed = 300",
    "10^\u22124 to 10^\u22121  (20 log-spaced values)",
    "3, 5, 7, 9, 11, 13, 15  [7 values]",
    "0, 1","TRUE (kernel density)","0.5, 1.0, 1.5",
    "50, 100, 150","2, 3, 4","0.05, 0.10","0, 0.1","Fixed = 0.80","Fixed = 1","Fixed = 0.80"
  ),
  Type = c(
    "Fixed",
    "Tuned","Fixed",
    "Tuned","Tuned",
    "Tuned","Fixed",
    "Tuned",
    "Tuned",
    "Tuned","Fixed","Tuned",
    "Tuned","Tuned","Tuned","Tuned","Fixed","Fixed","Fixed"
  ),
  stringsAsFactors = FALSE
)

model_fills_clf <- c(
  "LDA"="#e0f2fe",     "PLS-DA"="#dbeafe",   "SVM-RBF"="#dcfce7",
  "RF"="#ffedd5",      "CART"="#fef9c3",      "kNN"="#fce7f3",
  "Naive Bayes"="#f3e8ff", "XGBoost"="#fee2e2"
)
row_fills <- c(
  rep(model_fills_clf["LDA"],        1),
  rep(model_fills_clf["PLS-DA"],     2),
  rep(model_fills_clf["SVM-RBF"],    2),
  rep(model_fills_clf["RF"],         2),
  rep(model_fills_clf["CART"],       1),
  rep(model_fills_clf["kNN"],        1),
  rep(model_fills_clf["Naive Bayes"],3),
  rep(model_fills_clf["XGBoost"],    7)
)

tbl_hs_clf <- tableGrob(hyper_summary_clf, rows=NULL,
                        theme=ttheme_minimal(base_size=8.8,
                                             core    = list(fg_params=list(col=COL_DARK, hjust=0, x=.02),
                                                            bg_params=list(fill=row_fills, col=NA)),
                                             colhead = list(fg_params=list(fontface="bold", col="white", hjust=.5, x=.5),
                                                            bg_params=list(fill=COL_MID, col=NA))))

grid.draw(editGrob(tbl_hs_clf,
                   vp=viewport(x=.5, y=.46, width=.94, height=.65)))

# Leyenda de colores
leg_x <- .06
for(mn in names(model_fills_clf)){
  grid.rect(x=unit(leg_x,"npc"), y=unit(.065,"npc"),
            width=unit(.085,"npc"), height=unit(.024,"npc"),
            gp=gpar(fill=model_fills_clf[mn], col="gray60", lwd=.5))
  grid.text(mn, x=unit(leg_x+.045,"npc"), y=unit(.065,"npc"),
            just=c("left","center"),
            gp=gpar(fontsize=7.2, col=COL_DARK))
  leg_x <- leg_x + .118
}
draw_page_number(pg)
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("5. Boruta Variable Selection")

y_cur <- draw_body(paste0(
  "Boruta was run on the training set for each of the ",
  length(scatter_opts)*length(deriv_opts)," preprocessing combinations. ",
  "The algorithm wraps a Random Forest and iteratively removes features whose importance ",
  "falls below the maximum importance of randomly permuted shadow features (maxRuns=500). ",
  "Confirmed variables are retained; all others are discarded before model training. ",
  "The table below summarises, for each preprocessing, the number of variables confirmed ",
  "and their wavenumber identifiers."), y=.87)

if(length(boruta_vars_log) > 0){
  bvl_df <- bind_rows(lapply(boruta_vars_log, function(x)
    data.frame(Preprocessing = x$Preprocessing,
               N_selected    = x$N_selected,
               N_total       = if(exists("X_raw")) ncol(X_raw) else NA,
               Pct_selected  = if(exists("X_raw"))
                 paste0(round(x$N_selected/ncol(X_raw)*100,1),"%")
               else "N/A",
               stringsAsFactors=FALSE)))
  
  tbl_bv <- tableGrob(bvl_df, rows=NULL,
                      theme=ttheme_minimal(base_size=9,
                                           core   =list(fg_params=list(col=COL_DARK)),
                                           colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
  grid.draw(editGrob(tbl_bv,
                     vp=viewport(x=.5, y=max(y_cur-.18,.42), width=.80, height=.35)))
} else {
  grid.text("[boruta_vars_log not found — add logging to Boruta loop]",
            x=.5, y=y_cur-.08,
            gp=gpar(fontsize=9,col="red",fontface="italic"))
}
draw_page_number(pg)

# ── P boruta: variable lists (one page per preprocessing with variables) ──
if(length(boruta_vars_log) > 0){
  # Split into chunks of ~8 preprocessing combos per page
  chunks <- split(boruta_vars_log,
                  ceiling(seq_along(boruta_vars_log)/8))
  for(ck in seq_along(chunks)){
    grid.newpage(); pg <- pg + 1L
    draw_section_bar(paste0("5. Boruta — Variable Lists (page ",ck,")"))
    yy <- .88
    for(entry in chunks[[ck]]){
      # Header row
      grid.rect(x=unit(.05,"npc"),y=unit(yy,"npc"),
                width=unit(.90,"npc"),height=unit(.032,"npc"),
                just=c("left","top"),
                gp=gpar(fill="#eef4fb",col="#b0c8e0",lwd=.6))
      grid.text(paste0(entry$Preprocessing,
                       "  \u2014  ",entry$N_selected," variables confirmed"),
                x=unit(.07,"npc"),y=unit(yy-.016,"npc"),
                just=c("left","center"),
                gp=gpar(fontsize=9.5,fontface="bold",col=COL_MID))
      yy <- yy - .038
      # Variable list wrapped
      var_text <- if(nchar(entry$Variables)>0) entry$Variables else "(all variables used)"
      vlines   <- strwrap(var_text, width=115)
      for(vl in vlines){
        grid.text(vl,x=unit(.07,"npc"),y=unit(yy,"npc"),
                  just=c("left","top"),
                  gp=gpar(fontsize=7.5,col="#444444",fontfamily="mono"))
        yy <- yy - .022
      }
      yy <- yy - .012
      if(yy < .12) break
    }
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# SECTION 6 — BORUTA VARIABLES ON RAW SPECTRA
# ══════════════════════════════════════════════════════════════
if(exists("boruta_png_map_clf") && length(boruta_png_map_clf)>0){
  
  # Introductory page
  grid.newpage(); pg<-pg+1L
  draw_section_bar("6. Boruta Variable Selection — Highlighted on Raw Spectra")
  
  draw_body(paste0(
    "For each of the ",length(boruta_png_map_clf)," preprocessing combinations, the figure shows ",
    "the raw training spectra (post-outlier removal, n=",n_tr," samples, no mathematical ",
    "transformation applied) with the Boruta-confirmed variables marked as translucent red ",
    "vertical lines (\u03b1=0.35, lwd=0.9). Spectra are coloured by class."),
    y=.86)
  
  draw_body(paste0(
    "The Boruta selection was performed on the preprocessed version of the data for each ",
    "combination, but the selected wavenumber positions are always visualised on the original ",
    "raw spectrum. This ensures consistent chemical interpretation: the same physical absorption ",
    "bands appear at identical x-axis positions across all preprocessing combinations, enabling ",
    "direct cross-comparison of which spectral regions are discriminant for the target classes."),
    y=.74)
  
  # Recuadro interpretativo
  grid.rect(x=unit(.05,"npc"), y=unit(.26,"npc"),
            width=unit(.90,"npc"), height=unit(.28,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill="#f0f6ff", col=COL_ACCENT, lwd=1))
  grid.text("How to read these figures:",
            x=unit(.08,"npc"), y=unit(.52,"npc"),
            just=c("left","center"),
            gp=gpar(fontsize=10, fontface="bold", col=COL_MID))
  hints <- c(
    "\u25B6  Coloured lines = raw training spectra (blue = Class 1, red = Class 2)",
    "\u25B6  Red vertical lines = wavenumber positions confirmed by Boruta for that preprocessing",
    "\u25B6  Dense red marks = many variables selected; sparse = highly localised discriminant bands",
    "\u25B6  Regions marked across multiple preprocessings = chemically robust discriminant bands",
    "\u25B6  Each page header shows preprocessing name and number of confirmed variables"
  )
  for(ii in seq_along(hints))
    grid.text(hints[ii],
              x=unit(.09,"npc"), y=unit(.485-(ii-1)*.042,"npc"),
              just=c("left","center"),
              gp=gpar(fontsize=9, col="#333333"))
  draw_page_number(pg)
  
  # Overview without marks
  grid.newpage(); pg<-pg+1L
  draw_section_bar("6. Raw Spectra Overview — Classes (no Boruta marks)")
  if(file.exists(raw_spectra_clf_path)){
    grid.draw(editGrob(insert_png(raw_spectra_clf_path),
                       vp=viewport(x=.5, y=.48, width=.93, height=.82)))
    grid.text(paste0("Figure. Raw training spectra — ",n_tr," samples (post-outlier). ",
                     "Colour = class membership."),
              x=.5, y=.062, gp=gp_caption)
  } else {
    draw_body("[Raw spectra overview PNG not found]", y=.5)
  }
  draw_page_number(pg)
  
  # One page per preprocessing
  fig_n_clf <- 10L
  for(pp_nm in names(boruta_png_map_clf)){
    grid.newpage(); pg<-pg+1L; fig_n_clf<-fig_n_clf+1L
    n_sel <- if(!is.na(boruta_png_map_clf[[pp_nm]]$n_vars))
      boruta_png_map_clf[[pp_nm]]$n_vars else "N/A"
    
    draw_section_bar(paste0("6. Boruta on Raw Spectra — ",pp_nm,
                            "  [",n_sel," variables]"))
    
    grid.draw(editGrob(insert_png(boruta_png_map_clf[[pp_nm]]$raw_boruta),
                       vp=viewport(x=.5, y=.48, width=.93, height=.82)))
    
    grid.text(
      paste0("Figure ",fig_n_clf,". Raw training spectra (n=",n_tr,") with ",n_sel,
             " Boruta-confirmed wavenumbers marked in red (translucent). ",
             "Selection derived from preprocessed data: ",pp_nm,
             ". Colours = class membership."),
      x=.5, y=.062, gp=gp_caption)
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# SECTION 7 — RESULTS TABLES
# ══════════════════════════════════════════════════════════════
if(exists("results_df") && nrow(results_df) > 0){
  
  # Top 20 by AUC
  grid.newpage(); pg <- pg + 1L
  draw_section_bar("7. Performance Results — Top 20 Combinations (by AUC Test)")
  
  top20 <- results_df %>%
    arrange(desc(AUC_Test), desc(F1_Test)) %>% head(20) %>%
    mutate(Rank=row_number(),
           AUC_Test     =round(AUC_Test,    3),
           Accuracy_Test=round(Accuracy_Test,3),
           F1_Test      =round(F1_Test,     3),
           Kappa_Test   =round(Kappa_Test,  3),
           Sensitivity_Test=round(Sensitivity_Test,3),
           Specificity_Test=round(Specificity_Test,3)) %>%
    select(Rank,Model,Preprocessing,Boruta,N_vars,
           AUC_Test,Accuracy_Test,F1_Test,Kappa_Test,
           Sensitivity_Test,Specificity_Test)
  
  tbl20 <- tableGrob(top20, rows=NULL,
                     theme=ttheme_minimal(base_size=7.5,
                                          core   =list(fg_params=list(col=COL_DARK)),
                                          colhead=list(fg_params=list(fontface="bold",col=COL_MID,fontsize=8))))
  grid.draw(editGrob(tbl20,
                     vp=viewport(x=.5,y=.47,width=.96,height=.76)))
  grid.text("Sorted by AUC Test (descending), then F1 Test.",
            x=.5,y=.065,gp=gp_caption)
  draw_page_number(pg)
  
  # Summary by model
  grid.newpage(); pg <- pg + 1L
  draw_section_bar("7. Performance Results — Summary by Model")
  
  ms_tbl <- results_df %>% group_by(Model) %>%
    summarise(N        =n(),
              AUC_Mean =round(mean(AUC_Test,     na.rm=TRUE),3),
              AUC_Max  =round(max(AUC_Test,      na.rm=TRUE),3),
              Acc_Mean =round(mean(Accuracy_Test, na.rm=TRUE),3),
              Acc_Max  =round(max(Accuracy_Test,  na.rm=TRUE),3),
              F1_Mean  =round(mean(F1_Test,       na.rm=TRUE),3),
              Kappa_Mean=round(mean(Kappa_Test,   na.rm=TRUE),3),
              .groups="drop") %>% arrange(desc(AUC_Mean))
  
  tbl_ms <- tableGrob(ms_tbl, rows=NULL,
                      theme=ttheme_minimal(base_size=9.5,
                                           core   =list(fg_params=list(col=COL_DARK)),
                                           colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
  grid.draw(editGrob(tbl_ms,
                     vp=viewport(x=.5,y=.77,width=.88,height=.30)))
  hrule(.61)
  
  img_top10 <- insert_png("Results/Other/Top10_Combinations.png")
  grid.draw(editGrob(img_top10,
                     vp=viewport(x=.5,y=.355,width=.90,height=.44)))
  grid.text("Figure 5. Top 10 combinations ranked by AUC Test.",
            x=.5,y=.065,gp=gp_caption)
  draw_page_number(pg)
}

# ══════════════════════════════════════════════════════════════
# SECTION 8 — HEATMAPS
# ══════════════════════════════════════════════════════════════
hm_list <- list(
  list(f="Results/Heatmaps/Metrics_Model_Preprocessing.png",
       c="Figure 6. Accuracy / F1 / AUC by Model vs Preprocessing (mean across Boruta options)."),
  list(f="Results/Heatmaps/Metrics_Model_Boruta.png",
       c="Figure 7. Metrics by Model vs Preprocessing + Boruta configuration."),
  list(f="Results/Heatmaps/Metrics_by_Model.png",
       c="Figure 8. Average metrics aggregated by model."),
  list(f="Results/Heatmaps/Metrics_by_Preprocessing.png",
       c="Figure 9. Average metrics aggregated by preprocessing method."),
  list(f="Results/Heatmaps/AUC_Model_Preprocessing.png",
       c="Figure 10. AUC Test heatmap: Model vs Preprocessing. Red <0.80; Yellow 0.80-0.95; Green \u22650.95."),
  list(f="Results/Other/Model_Comparison_AUC.png",
       c="Figure 11. AUC Test by model and preprocessing. Points = individual combinations; bar = median."),
  list(f="Results/Other/Model_Comparison_Accuracy.png",
       c="Figure 12. Accuracy Test by model and preprocessing."),
  list(f="Results/Other/PCA_Metrics.png",
       c="Figure 13. PCA of all performance metrics (Train + Test). Point size \u221d AUC Test.")
)
for(hm in hm_list){
  grid.newpage(); pg <- pg + 1L
  draw_section_bar("8. Heatmaps — Model vs Preprocessing")
  grid.draw(editGrob(insert_png(hm$f),
                     vp=viewport(x=.5,y=.495,width=.92,height=.82)))
  grid.text(hm$c,x=.5,y=.062,gp=gp_caption)
  draw_page_number(pg)
}

# ══════════════════════════════════════════════════════════════
# SECTION 9 — CONFUSION MATRICES  (2 per page: Train + Test)
# ══════════════════════════════════════════════════════════════
clean_cap <- function(fname){
  b <- tools::file_path_sans_ext(basename(fname))
  b <- gsub("_"," ",b); trimws(gsub("  +"," ",b))
}

for(m_name in c(names(models),"XGB")){
  cm_dir <- paste0("Results/ConfusionMatrices/",m_name)
  if(!dir.exists(cm_dir)) next
  all_cm <- list.files(cm_dir, pattern="\\.png$", full.names=TRUE)
  if(length(all_cm)==0) next
  
  # Header page for this model
  grid.newpage(); pg <- pg + 1L
  draw_section_bar(paste0("9. Confusion Matrices — ",m_name))
  
  if(exists("results_df") && nrow(results_df)>0){
    best_m <- results_df %>% filter(Model==m_name) %>%
      arrange(desc(AUC_Test), desc(Accuracy_Test), desc(Kappa_Test), N_vars) %>% head(5) %>%
      mutate(AUC_Test     =round(AUC_Test,    3),
             Accuracy_Test=round(Accuracy_Test,3),
             F1_Test      =round(F1_Test,     3),
             Kappa_Test   =round(Kappa_Test,  3)) %>%
      select(Preprocessing,Boruta,N_vars,AUC_Test,Accuracy_Test,F1_Test,Kappa_Test)
    tbl_b <- tableGrob(best_m, rows=NULL,
                       theme=ttheme_minimal(base_size=9,
                                            core   =list(fg_params=list(col=COL_DARK)),
                                            colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
    grid.draw(editGrob(tbl_b,
                       vp=viewport(x=.5,y=.55,width=.80,height=.28)))
    grid.text(paste0("Top 5 combinations for ",m_name," sorted by AUC Test."),
              x=.5,y=.38,gp=gp_caption)
  } else {
    draw_body(paste0("Confusion matrices for all preprocessing × Boruta combinations ",
                     "evaluated with ",m_name,". Train (left) and Test (right) shown side by side."),
              y=.75)
  }
  draw_page_number(pg)
  
  # Pair Train + Test side by side on each page
  train_files <- sort(all_cm[grepl("Train",all_cm)])
  test_files  <- sort(all_cm[grepl("Test", all_cm)])
  n_pairs <- max(length(train_files), length(test_files))
  
  for(pi in seq_len(n_pairs)){
    grid.newpage(); pg <- pg + 1L
    grid.text(paste0(m_name," — Page ",pi," of ",n_pairs),
              x=.5,y=.974,
              gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
    grid.lines(x=c(.04,.96),y=c(.955,.955),gp=gpar(col="gray75",lwd=.6))
    
    if(pi <= length(train_files)){
      grid.draw(editGrob(insert_png(train_files[pi]),
                         vp=viewport(x=.25,y=.50,width=.46,height=.82)))
      grid.text(paste0("TRAIN — ",clean_cap(train_files[pi])),
                x=.25,y=.075,gp=gp_caption)
    }
    if(pi <= length(test_files)){
      grid.draw(editGrob(insert_png(test_files[pi]),
                         vp=viewport(x=.75,y=.50,width=.46,height=.82)))
      grid.text(paste0("TEST — ",clean_cap(test_files[pi])),
                x=.75,y=.075,gp=gp_caption)
    }
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# SECTION 10 — ROC CURVES  (2 per page)
# ══════════════════════════════════════════════════════════════
for(m_name in c(names(models),"XGB")){
  roc_dir <- paste0("Results/ROC/",m_name)
  if(!dir.exists(roc_dir)) next
  roc_files <- list.files(roc_dir, pattern="\\.png$", full.names=TRUE)
  if(length(roc_files)==0) next
  
  # header page
  grid.newpage(); pg <- pg + 1L
  draw_section_bar(paste0("10. ROC Curves — ",m_name))
  draw_body(paste0(
    "ROC curves for all preprocessing \u00d7 Boruta combinations evaluated with ",
    m_name,". Each panel shows the test-set ROC curve with AUC. ",
    "Positive class: ",if(exists("POS_CLASS")) POS_CLASS else "N/A","."),
    y=.78)
  draw_page_number(pg)
  
  # 2 ROC curves per page
  n_roc <- length(roc_files)
  n_pg_r <- ceiling(n_roc/2)
  for(rp in seq_len(n_pg_r)){
    grid.newpage(); pg <- pg + 1L
    grid.text(paste0(m_name," ROC — Page ",rp," of ",n_pg_r),
              x=.5,y=.974,
              gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
    grid.lines(x=c(.04,.96),y=c(.955,.955),gp=gpar(col="gray75",lwd=.6))
    
    idx1 <- (rp-1)*2+1; idx2 <- min(rp*2, n_roc)
    grid.draw(editGrob(insert_png(roc_files[idx1]),
                       vp=viewport(x=.27,y=.50,width=.50,height=.84)))
    grid.text(clean_cap(roc_files[idx1]),x=.27,y=.065,gp=gp_caption)
    
    if(idx2 > idx1){
      grid.draw(editGrob(insert_png(roc_files[idx2]),
                         vp=viewport(x=.77,y=.50,width=.50,height=.84)))
      grid.text(clean_cap(roc_files[idx2]),x=.77,y=.065,gp=gp_caption)
    }
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# CLOSE PDF
# ══════════════════════════════════════════════════════════════
dev.off()
cat("\n",rep("=",60),"\n",sep="")
cat("  CLASSIFICATION PDF REPORT GENERATED\n")
cat("  Path:  ",pdf_path,"\n",sep="")
cat("  Pages: ",pg,"\n",sep="")
cat(rep("=",60),"\n\n",sep="")

# ══════════════════════════════════════════════════════════════
# Rename the results folder using the name of the analyzed sheet
# (e.g. "Results" -> "Results_HIERRO"), so that the results of
# another analyte are not overwritten in the next run.
# ══════════════════════════════════════════════════════════════
results_root <- paste0("Results_", sheet_name)
if (dir.exists(results_root)) {
  cat(">>> '",results_root,"' already existed (previous run for the same analyte) - it will be overwritten.\n",sep="")
  unlink(results_root, recursive = TRUE)
}
if (dir.exists("Results")) {
  file.rename("Results", results_root)
  cat(">>> Results folder renamed to: ",results_root,"\n",sep="")
}

# ══════════════════════════════════════════════════════════════
# Append the analyte name to ALL generated files
# (ej. "Boruta_Paired_Ttest.xlsx" -> "Boruta_Paired_Ttest_POTASSIUM.xlsx").
# Automatic: the suffix is taken from the name of the analyzed sheet (sheet_name).
# Done at the end, once all files are closed, so as not to interfere
# with the rest of the pipeline (plots and the PDF are generated with the
# original names and only renamed afterwards).
# ══════════════════════════════════════════════════════════════
if (dir.exists(results_root)) {
  tag   <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
  files <- list.files(results_root, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  base  <- basename(files)
  todo  <- grepl("\\.[^.]+$", base) & !grepl(paste0("_", tag, "\\.[^.]+$"), base)
  new   <- file.path(dirname(files), sub("(\\.[^.]+)$", paste0("_", tag, "\\1"), base))
  ok    <- rep(FALSE, length(files))
  if (any(todo)) ok[todo] <- file.rename(files[todo], new[todo])
  cat(">>> Files renamed with the suffix '_", tag, "': ", sum(ok), " of ", sum(todo), "\n", sep = "")
  if (sum(ok) < sum(todo))
    cat("    (files that could not be renamed are usually open in another program, e.g. Excel or the PDF viewer)\n")
}
