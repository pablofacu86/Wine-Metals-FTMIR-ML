# ======================================================================
# PIPELINE SCAN REGRESSION
# Pipeline for the quantitative estimation of potassium, magnesium and calcium
# in wine from FT-MIR spectra. See the README for a complete description of the
# workflow, the input data and the generated outputs.
# ======================================================================


# ======================================
# PREPROCESSING SETTINGS
# ======================================
# SCATTER_FIRST: order of operations, always applied to the COMPLETE spectrum.
#   TRUE  -> scatter correction (SNV/MSC) -> Savitzky-Golay (smoothing or derivative)
#            SNV and MSC were designed for underived absorbance spectra.
#   FALSE -> Savitzky-Golay -> SNV/MSC (alternative order).
# SMOOTH_ALWAYS: if TRUE, combinations WITHOUT a derivative are also smoothed with
#   Savitzky-Golay (derivative order 0, polynomial 3, 11-point window). If FALSE,
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
cat("FT-MIR REGRESSION PIPELINE (Scan Regression)\n")
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
library(pls)
library(e1071)
library(randomForest)
library(openxlsx)
library(prospectr)
library(Boruta)
library(xgboost)
library(ggplot2)
library(glmnet)
library(ggrepel)
library(ggnewscale)
library(scales)

# Resolve name conflict for the margin function
margin <- ggplot2::margin

# ======================================
# GLOBAL SETTINGS
# ======================================
CONC_UNITS <- "mg/L"   # concentration units of the analyte

# ======================================
# 2. Directories
# ======================================
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

dir.create("Results",               showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Excel",         showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Spectra",       showWarnings = FALSE, recursive = TRUE)
dir.create("Results/ScatterPlots",  showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Heatmaps",      showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Other",         showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Outliers",      showWarnings = FALSE, recursive = TRUE)

cat(">>> Output folders created.\n\n")

# ======================================
# 3. Data loading
# ======================================
cat(">>> Loading data...\n")

data_path  <- "data/FINAL_DATA_SET.xlsx"      # path to the data file (see the README)
sheet_name <- "MAGNESIUM"   # Analyte to run: "POTASSIUM", "MAGNESIUM" or "CALCIUM"

# Locates the data file regardless of the working directory and of whether the
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

nir <- read.xlsx(data_path, sheet = sheet_in_file)
nir <- data.frame(lapply(nir, function(x)
  if (is.character(x)) as.factor(x) else x))

sample_ids  <- as.character(nir[, 1])
y_raw       <- nir[, 2]
X_raw       <- nir[, -c(1, 2)]                     # COMPLETE, continuous spectrum

# --------------------------------------------------------------------
# Helper: parses wavenumbers from column names such as "X960.45"
# or "X960,45" (decimal comma, regional Excel format). Without this
# replacement, as.numeric("960,45") silently returns NA and
# breaks patz_idx, the post-preprocessing cropping and the plots with
# Boruta lines. It is used throughout the pipeline instead of a bare
# as.numeric(gsub("X","",...)).
# --------------------------------------------------------------------
parse_wn <- function(x) as.numeric(gsub(",", ".", gsub("X", "", x), fixed = TRUE))

# --------------------------------------------------------------------
# Helpers: format the optimized hyperparameters of each winning model
# (caret$bestTune, or the grid selected by manual CV in XGB)
# as a readable string, so they can be reported together with the
# performance metrics of each combination.
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
cat("   Response range:", round(min(y_raw), 2), "-", round(max(y_raw), 2), "\n\n")

# --------------------------------------------------------------------
# Spectral windows from Patz et al. (2004), retained for modeling.
# The spectrum is cropped AFTER smoothing/derivation, to avoid artifacts
# at the junction of non-contiguous windows.
# --------------------------------------------------------------------
PATZ_WINDOWS <- list(c(965, 1582), c(1698, 2006), c(2701, 2971))
in_patz_windows <- function(w) {
  Reduce(`|`, lapply(PATZ_WINDOWS, function(rng) w >= rng[1] & w <= rng[2]))
}
patz_idx <- in_patz_windows(wavelengths)
cat("   Spectral variables within Patz windows:", sum(patz_idx), "\n\n")

# ======================================
# 4. Outlier detection and removal
#    (on the cropped Patz variables)
# ======================================
cat(">>> Outlier detection (Hotelling T2 AND Q residuals)...\n")

X_raw_patz <- X_raw[, patz_idx, drop = FALSE]
X_scaled <- scale(X_raw_patz)
pca_exp  <- prcomp(X_scaled, center = FALSE, scale. = FALSE)
var_cum  <- cumsum(pca_exp$sdev^2) / sum(pca_exp$sdev^2)
n_pcs    <- max(3, which(var_cum >= 0.95)[1])
cat("   PCs used:", n_pcs, sprintf("(%.1f%% variance)\n", var_cum[n_pcs] * 100))

scores   <- pca_exp$x[, 1:n_pcs, drop = FALSE]
var_exp2 <- round(pca_exp$sdev^2 / sum(pca_exp$sdev^2) * 100, 1)

lambda_inv <- diag(1 / pca_exp$sdev[1:n_pcs]^2)
T2         <- rowSums((scores %*% lambda_inv) * scores)
T2_lim     <- qchisq(0.99, df = n_pcs)
flag_T2    <- T2 > T2_lim

X_rec  <- scores %*% t(pca_exp$rotation[, 1:n_pcs, drop = FALSE])
Q      <- rowSums((X_scaled - X_rec)^2)
Q_lim  <- mean(Q) + 3 * sd(Q)
flag_Q <- Q > Q_lim

flag_outlier <- flag_T2 & flag_Q   # AND: only "severe" outliers that violate BOTH criteria at once
n_outliers   <- sum(flag_outlier)

cat("   T2 flagged:", sum(flag_T2), "| Q flagged:", sum(flag_Q),
    "| Either (removed):", n_outliers, "of", nrow(X_raw), "samples\n")

if (n_outliers > 0) {
  oidx <- which(flag_outlier)
  write.xlsx(data.frame(
    Sample_ID=sample_ids[oidx], Row_Index=oidx,
    Y_Value=as.numeric(y_raw)[oidx],
    T2_score=round(T2[oidx],3),   T2_limit=round(T2_lim,3),   T2_flagged=flag_T2[oidx],
    Q_residual=round(Q[oidx],4),  Q_limit=round(Q_lim,4),     Q_flagged=flag_Q[oidx],
    Criterion="T2 > chi2(0.99) AND Q > mean+3SD", stringsAsFactors=FALSE),
    file="Results/Outliers/Outliers_Removed.xlsx", overwrite=TRUE)
  cat("   Removed:", paste(sample_ids[oidx], collapse=", "), "\n")
} else {
  write.xlsx(data.frame(Message="No outliers detected - all samples retained.",
    T2_limit=round(T2_lim,3), Q_limit=round(Q_lim,4),
    Criterion="T2 AND Q", stringsAsFactors=FALSE),
    file="Results/Outliers/Outliers_Removed.xlsx", overwrite=TRUE)
  cat("   No outliers detected.\n")
}

# ── Raw spectra ─────────────────────────────────────────────────
cat("   Plotting raw spectra...\n")
raw_long <- data.frame(
  Wavenumber=rep(wavelengths, times=nrow(X_raw)),
  Absorbance=as.vector(t(as.matrix(X_raw))),
  Sample_ID =rep(sample_ids, each=length(wavelengths)),
  Group     =rep(ifelse(flag_outlier,"Outlier","Normal"), each=length(wavelengths)),
  stringsAsFactors=FALSE)

p_raw <- ggplot() +
  geom_line(data=subset(raw_long,Group=="Normal"),
            aes(x=Wavenumber,y=Absorbance,group=Sample_ID),
            color="#2171B5",alpha=0.35,linewidth=0.4) +
  geom_line(data=subset(raw_long,Group=="Outlier"),
            aes(x=Wavenumber,y=Absorbance,group=Sample_ID),
            color="#D7191C",alpha=0.85,linewidth=0.7) +
  labs(title="Raw FTMIR Spectra",
       subtitle=paste0("Blue: normal (n=",sum(!flag_outlier),")  |  Red: outliers (n=",n_outliers,")"),
       x=expression("Wavenumber (cm"^{-1}*")"), y="Absorbance") +
  theme_minimal(base_size=12) +
  theme(plot.title=element_text(face="bold",hjust=0.5,size=14),
        plot.subtitle=element_text(hjust=0.5,size=9,color="gray40"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="gray70",fill=NA))
if (n_outliers>0) {
  olbl <- subset(raw_long,Group=="Outlier") %>%
    group_by(Sample_ID) %>% slice_max(Absorbance,n=1) %>% ungroup()
  p_raw <- p_raw +
    geom_label_repel(data=olbl, aes(x=Wavenumber,y=Absorbance,label=Sample_ID),
                     color="#D7191C",fill="white",fontface="bold",
                     size=3,box.padding=0.4,max.overlaps=20,
                     segment.color="#D7191C",segment.linewidth=0.4)
}
ggsave("Results/Outliers/Raw_Spectra_Outliers.png",p_raw,width=12,height=6,dpi=300,bg="white")
cat("   OK Raw spectra saved\n")

# ── PCA plots ────────────────────────────────────────────────────────
pca_df <- data.frame(PC1=pca_exp$x[,1],PC2=pca_exp$x[,2],PC3=pca_exp$x[,3],
                     Sample_ID=sample_ids,Y_value=as.numeric(y_raw),
                     Outlier=ifelse(flag_outlier,"Outlier","Normal"),stringsAsFactors=FALSE)

make_pca_plot <- function(df,xc,yc,xlab,ylab,title,fname) {
  dn<-subset(df,Outlier=="Normal"); do<-subset(df,Outlier=="Outlier")
  p<-ggplot(df,aes_string(x=xc,y=yc))+
    geom_point(data=dn,aes(color=Y_value),size=3.5,alpha=0.85)+
    scale_color_viridis_c(name="Y value",option="D")+
    geom_hline(yintercept=0,linetype="dashed",color="gray60",linewidth=0.4)+
    geom_vline(xintercept=0,linetype="dashed",color="gray60",linewidth=0.4)+
    labs(title=title,x=xlab,y=ylab)+
    theme_minimal(base_size=12)+
    theme(plot.title=element_text(face="bold",hjust=0.5,size=13),
          legend.position="right",panel.border=element_rect(color="gray70",fill=NA))
  if(nrow(do)>0)
    p<-p+geom_point(data=do,color="red",fill="red",shape=23,size=4.5,stroke=1.2)+
      annotate("text",x=min(df[[xc]],na.rm=TRUE),y=max(df[[yc]],na.rm=TRUE),
               label=paste0("Outliers removed: ",nrow(do)),
               hjust=0,vjust=1,color="red",fontface="bold",size=3.8)
  ggsave(fname,p,width=11,height=8,dpi=300,bg="white")
  cat("   OK PCA saved:",basename(fname),"\n"); invisible(p)
}
make_pca_plot(pca_df,"PC1","PC2",
              paste0("PC1 (",var_exp2[1],"%)"),paste0("PC2 (",var_exp2[2],"%)"),
              "Exploratory PCA - PC1 vs PC2","Results/Outliers/PCA_PC1_vs_PC2.png")
make_pca_plot(pca_df,"PC2","PC3",
              paste0("PC2 (",var_exp2[2],"%)"),paste0("PC3 (",var_exp2[3],"%)"),
              "Exploratory PCA - PC2 vs PC3","Results/Outliers/PCA_PC2_vs_PC3.png")

infl_df<-data.frame(T2=T2,Q=Q,Sample_ID=sample_ids,
                     Outlier=ifelse(flag_outlier,"Outlier","Normal"),stringsAsFactors=FALSE)
p_infl<-ggplot(infl_df,aes(x=T2,y=Q,color=Outlier))+
  geom_point(size=4,alpha=0.85)+
  geom_text_repel(data=subset(infl_df,Outlier=="Normal"),aes(label=Sample_ID),
                  size=3.4,color="gray40",box.padding=0.2,max.overlaps=20)+
  geom_vline(xintercept=T2_lim,linetype="dashed",color="firebrick",linewidth=0.8)+
  geom_hline(yintercept=Q_lim,linetype="dashed",color="steelblue",linewidth=0.8)+
  scale_color_manual(values=c("Normal"="steelblue","Outlier"="red"))+
  labs(title="Influence Plot: Hotelling T2 vs Q",
       subtitle="Outliers: T2 > chi2(0.99) AND Q > mean+3SD",
       x="Hotelling T2",y="Q Residuals",color="")+
  theme_minimal(base_size=22)+
  theme(plot.title=element_text(face="bold",hjust=0.5,size=24),
        plot.subtitle=element_text(hjust=0.5,size=17,color="gray40"),
        axis.title=element_text(size=20),axis.text=element_text(size=17),
        legend.text=element_text(size=17),legend.title=element_text(size=18),
        legend.position="top",panel.border=element_rect(color="gray70",fill=NA))
if(n_outliers>0)
  p_infl<-p_infl+
    geom_label_repel(data=subset(infl_df,Outlier=="Outlier"),aes(label=Sample_ID),
                     color="red",fill="white",fontface="bold",
                     size=4.5,box.padding=0.4,max.overlaps=30)
ggsave("Results/Outliers/Influence_Plot_T2_vs_Q.png",p_infl,width=11,height=8,dpi=300,bg="white")
cat("   OK Influence plot saved\n")

keep_idx<-!flag_outlier; X<-X_raw[keep_idx,,drop=FALSE]
y<-y_raw[keep_idx]; sample_ids<-sample_ids[keep_idx]
cat(sprintf("\n   After removal: %d samples (removed %d)\n\n",sum(keep_idx),n_outliers))

# ======================================
# 5. Train / Test Split (stratified random)
# ======================================
# Stratified random sampling (createDataPartition) with a fixed seed.
# Distance-based algorithms (Kennard-Stone, Duplex) are avoided, since they
# concentrate the extreme samples in one of the two sets.
cat(">>> Splitting data (70/30, stratified random sampling over 5 concentration bins)...\n")
set.seed(1234)
bins<-cut(y,breaks=5); idx<-createDataPartition(bins,p=0.7,list=FALSE)

X_train<-X[idx,]; X_test<-X[-idx,]; y_train<-y[idx]; y_test<-y[-idx]
cat("   Train:",nrow(X_train),"| Test:",nrow(X_test),"\n\n")

# ======================================
# 6. Analysis functions
# ======================================
apply_preprocessing<-function(Xtr,Xte,ytr,yte,scatter,deriv){
  tryCatch({
    Xtr<-as.matrix(Xtr); Xte<-as.matrix(Xte)
    co<-colnames(Xtr); wo<-parse_wn(co)

    # --- Order of operations (SCATTER_FIRST / SMOOTH_ALWAYS switches) ---
    # Everything is computed on the COMPLETE, continuous spectrum (545 variables, no
    # gaps), so the 11-point Savitzky-Golay moving window does not mix signal
    # from non-contiguous regions. Cropping to the Patz windows is done afterwards.
    keep_names <- function(Xnew, cn) {
      # Savitzky-Golay may trim (w-1)/2 columns at each edge; the remaining
      # columns are the central ones, so they are named by position.
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
        ref <- colMeans(Xa, na.rm = TRUE)   # reference: mean of the TRAINING set
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

    # --- Crop to the Patz windows AFTER smoothing/derivative/scatter correction ---
    co2<-colnames(Xtr); wo2<-suppressWarnings(parse_wn(co2))
    if(length(wo2)==0 || all(is.na(wo2))) wo2<-wo   # fallback in case the names were lost
    keep_patz<-in_patz_windows(wo2)
    Xtr<-Xtr[,keep_patz,drop=FALSE]; Xte<-Xte[,keep_patz,drop=FALSE]
    wo2<-wo2[keep_patz]

    mu<-colMeans(Xtr,na.rm=TRUE)
    Xtr_plot<-Xtr; Xte_plot<-Xte   # <- UNcentered version, only for diagnostic plots
    Xtr<-sweep(Xtr,2,mu,"-");Xte<-sweep(Xte,2,mu,"-")
    ktr<-complete.cases(Xtr);kte<-complete.cases(Xte)
    vv<-apply(Xtr[ktr,,drop=FALSE],2,var,na.rm=TRUE);kc<-!is.na(vv)&vv>1e-10
    fn<-colnames(Xtr)[kc];fw<-suppressWarnings(parse_wn(fn))
    if(any(is.na(fw)))fw<-wo2[kc]
    list(Xtr=Xtr[ktr,kc,drop=FALSE],Xte=Xte[kte,kc,drop=FALSE],
         Xtr_plot=Xtr_plot[ktr,kc,drop=FALSE],
         ytr=ytr[ktr],yte=yte[kte],wavelengths=fw,success=TRUE)
  },error=function(e)list(success=FALSE,error=as.character(e)))
}

plot_spectra<-function(Xm,wl,title,sel_vars=NULL,pp_name){
  nc<-ncol(Xm)
  if(length(wl)>nc)wl<-wl[1:nc] else if(length(wl)<nc)wl<-1:nc
  if(nrow(Xm)==0||nc==0)return(NULL)
  df<-data.frame(wavelength=rep(wl,each=nrow(Xm)),
                  absorbance=as.vector(Xm),sample=rep(1:nrow(Xm),times=length(wl)))
  df<-df[complete.cases(df),];if(nrow(df)==0)return(NULL)
  p<-ggplot(df,aes(x=wavelength,y=absorbance,group=sample))+
    geom_line(alpha=0.3,color="gray40")+
    labs(title=title,x=expression("Wavenumber (cm"^{-1}*")"),y="Absorbance")+
    theme_minimal(base_size=15)+
    theme(plot.title=element_text(face="bold",hjust=0.5,size=16),
          axis.title=element_text(size=14),axis.text=element_text(size=12),
          panel.grid.minor=element_blank())
  if(!is.null(sel_vars)){
    sw<-parse_wn(sel_vars);sw<-sw[sw%in%wl]
    if(length(sw)>0)
      p<-p+geom_vline(xintercept=sw,color="red",alpha=0.6,linetype="dashed")+
        annotate("text",x=min(wl),y=max(df$absorbance,na.rm=TRUE),
                 label=paste0(length(sw)," variables selected"),
                 hjust=0,vjust=1,color="red",size=5,fontface="bold")
  }
  fn<-paste0("Results/Spectra/",gsub(" ","_",pp_name),
              ifelse(is.null(sel_vars),"","_Boruta"),".png")
  tryCatch(ggsave(fn,p,width=10,height=6,dpi=300,bg="white"),error=function(e)NULL)
  invisible(p)
}

calc_metrics<-function(model,tr,te){
  tryCatch({
    pt<-as.numeric(predict(model,tr));pe<-as.numeric(predict(model,te))
    yt<-tr$Grupo;ye<-te$Grupo
    list(R2_Train=caret::R2(pt,yt),RMSE_Train=caret::RMSE(pt,yt),RPD_Train=sd(yt)/caret::RMSE(pt,yt),
         R2_Test =caret::R2(pe,ye),RMSE_Test =caret::RMSE(pe,ye),RPD_Test =sd(ye)/caret::RMSE(pe,ye),
         success=TRUE)
  },error=function(e)list(success=FALSE,error=as.character(e)))
}

plot_predictions<-function(model,train_data,test_data,model_name,pp_name,boruta_status){
  tryCatch({
    pt<-as.numeric(predict(model,train_data));pe<-as.numeric(predict(model,test_data))
    rt<-train_data$Grupo;re<-test_data$Grupo
    r2t<-caret::R2(pt,rt);rmset<-caret::RMSE(pt,rt);rpdt<-sd(rt)/rmset
    r2e<-caret::R2(pe,re);rmsee<-caret::RMSE(pe,re);rpde<-sd(re)/rmsee
    pd<-data.frame(Real=c(rt,re),Predicted=c(pt,pe),
                   Type=c(rep("Training",length(rt)),rep("Test",length(re))))
    mn<-min(c(pd$Real,pd$Predicted),na.rm=TRUE);mx<-max(c(pd$Real,pd$Predicted),na.rm=TRUE)
    p<-ggplot(pd,aes(x=Real,y=Predicted,color=Type))+
      geom_point(size=3.5,alpha=0.7)+
      geom_abline(intercept=0,slope=1,linetype="dashed",color="black",linewidth=1)+
      scale_color_manual(values=c("Training"="steelblue","Test"="firebrick"))+
      labs(title=paste0(model_name," - ",pp_name),subtitle=paste0("Boruta: ",boruta_status),
           x="Actual Concentration",y="Predicted Concentration",color="Dataset")+
      annotate("text",x=mn,y=mx,
               label=paste0("TRAINING\nR2=",round(r2t,3),"\nRMSE=",round(rmset,3),"\nRPD=",round(rpdt,2)),
               hjust=0,vjust=1,size=4.3,color="steelblue",fontface="bold")+
      annotate("text",x=mx,y=mn,
               label=paste0("TEST\nR2=",round(r2e,3),"\nRMSE=",round(rmsee,3),"\nRPD=",round(rpde,2)),
               hjust=1,vjust=0,size=4.3,color="firebrick",fontface="bold")+
      theme_minimal(base_size=16)+
      theme(plot.title=element_text(face="bold",hjust=0.5,size=18),
            plot.subtitle=element_text(hjust=0.5,size=13),
            axis.title=element_text(size=15),axis.text=element_text(size=13),
            legend.text=element_text(size=13),legend.title=element_text(size=14),
            legend.position="top")+coord_fixed()
    dir.create(paste0("Results/ScatterPlots/",model_name),showWarnings=FALSE,recursive=TRUE)
    fn<-paste0("Results/ScatterPlots/",model_name,"/",
               gsub(" ","_",pp_name),"_",boruta_status,".png")
    ggsave(fn,p,width=8,height=8,dpi=300,bg="white");return(TRUE)
  },error=function(e){cat("   ERROR scatter:",e$message,"\n");return(FALSE)})
}

# ======================================================
# HEATMAP TYPE 1: facet_grid (hm1 and hm2)
#
# pivot_longer is used to create
# a long table with columns (x_col, Model, Metric, Value)
# and Metric is mapped to the Y axis within each facet.
# facet_grid(Model ~ .) groups the 3 metric rows
# under each model panel. The legends are
# native ggplot guide_colorbar objects, combined with ggnewscale.
# ======================================================
make_heatmap_facet <- function(df_wide, x_var, y_var="Model",
                                title, subtitle, filename,
                                width=18, height=11) {

  rmse_vals <- df_wide$RMSE_Test
  q33  <- quantile(rmse_vals, 1/3, na.rm=TRUE)
  q66  <- quantile(rmse_vals, 2/3, na.rm=TRUE)
  rmin <- min(rmse_vals, na.rm=TRUE)
  rmax <- max(rmse_vals, na.rm=TRUE)
  rpd_max <- max(df_wide$RPD_Test, na.rm=TRUE)

  lbl_rmse <- paste0("RMSE Test (", CONC_UNITS, ")")
  lbl_rpd  <- "RPD Test"
  lbl_r2   <- "R\u00b2 Test"

  # Long table: one row per (x_col, Model, Metric)
  # Y axis = Metric, facet = Model (or y_var)
  df_long <- df_wide %>%
    rename(x_col = !!sym(x_var), Model = !!sym(y_var)) %>%
    pivot_longer(cols = c(RMSE_Test, RPD_Test, R2_Test),
                 names_to  = "metric_id",
                 values_to = "Value") %>%
    mutate(Metric = factor(
      case_when(metric_id == "RMSE_Test" ~ lbl_rmse,
                metric_id == "RPD_Test"  ~ lbl_rpd,
                metric_id == "R2_Test"   ~ lbl_r2),
      levels = c(lbl_rmse, lbl_rpd, lbl_r2)   # order: RMSE on top, R2 at the bottom
    ))

  d_rmse <- filter(df_long, metric_id == "RMSE_Test")
  d_rpd  <- filter(df_long, metric_id == "RPD_Test")
  d_r2   <- filter(df_long, metric_id == "R2_Test")

  p <- ggplot() +

    # ── RMSE (green=low, red=high) ──────────────────────────────
    geom_tile(data=d_rmse, aes(x=x_col, y=Metric, fill=Value),
              color="white", linewidth=0.5) +
    geom_text(data=d_rmse, aes(x=x_col, y=Metric, label=round(Value,2)),
              color="black", fontface="bold", size=2.7) +
    scale_fill_gradientn(
      name    = paste0("RMSE\n(", CONC_UNITS, ")"),
      colours = c("darkgreen","yellowgreen","yellow","orange","red"),
      values  = rescale(c(rmin, q33, (q33+q66)/2, q66, rmax)),
      limits  = c(rmin, rmax),
      guide   = guide_colorbar(
        barheight=unit(2.5,"cm"), barwidth=unit(0.4,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=8,face="bold",color="black"),
        label.theme=element_text(size=7,color="black"),
        order=3)
    ) +
    new_scale_fill() +

    # ── RPD (red=low, green=high) ───────────────────────────────
    geom_tile(data=d_rpd, aes(x=x_col, y=Metric, fill=Value),
              color="white", linewidth=0.5) +
    geom_text(data=d_rpd, aes(x=x_col, y=Metric, label=round(Value,2)),
              color="black", fontface="bold", size=2.7) +
    scale_fill_gradientn(
      name    = "RPD Test",
      colours = c("red","orange","yellow","yellowgreen","darkgreen"),
      values  = rescale(c(0, 1.0, 1.50, 2.50, rpd_max)),
      limits  = c(0, rpd_max),
      guide   = guide_colorbar(
        barheight=unit(2.5,"cm"), barwidth=unit(0.4,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=8,face="bold",color="black"),
        label.theme=element_text(size=7,color="black"),
        order=2)
    ) +
    new_scale_fill() +

    # ── R² (red=low, green=high) ────────────────────────────────
    geom_tile(data=d_r2, aes(x=x_col, y=Metric, fill=Value),
              color="white", linewidth=0.5) +
    geom_text(data=d_r2, aes(x=x_col, y=Metric, label=round(Value,3)),
              color="black", fontface="bold", size=2.7) +
    scale_fill_gradientn(
      name    = "R\u00b2 Test",
      colours = c("red","orange","yellow","yellowgreen","darkgreen"),
      values  = rescale(c(0, 0.60, 0.80, 0.90, 1.0)),
      limits  = c(0, 1),
      guide   = guide_colorbar(
        barheight=unit(2.5,"cm"), barwidth=unit(0.4,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=8,face="bold",color="black"),
        label.theme=element_text(size=7,color="black"),
        order=1)
    ) +

    # facet_grid: Model in rows, the 3 metrics are the Y axis of each panel
    facet_grid(rows   = vars(Model),
               scales = "free_y",
               space  = "free_y",
               switch = "y") +

    labs(title=title, subtitle=subtitle, x=x_var, y=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_minimal(base_size=9.5) +
    theme(
      plot.title    = element_text(face="bold", hjust=0.5, size=13),
      plot.subtitle = element_text(hjust=0.5, size=7.5, color="gray40",
                                   margin=margin(b=6)),
      strip.placement   = "outside",
      strip.text.y.left = element_text(angle=0, face="bold", size=9,
                                        hjust=0.5, vjust=0.5,
                                        margin=margin(r=5,l=5)),
      strip.background  = element_rect(fill="gray88", color="black", linewidth=0.9),
      panel.border      = element_rect(color="black", fill=NA, linewidth=0.9),
      panel.spacing.y   = unit(2,"pt"),
      axis.text.x  = element_text(angle=45, hjust=1, size=7.5, color="black"),
      axis.text.y  = element_text(size=8, color="black"),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(8,5,8,5),
      legend.position   = "right",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.5,"cm"),
      legend.key.size   = unit(0.4,"cm")
    )

  ggsave(filename, p, width=width, height=height, dpi=300,bg="white")
  cat("   OK Heatmap saved:", basename(filename), "\n")
  invisible(p)
}

# ======================================================
# HEATMAP TYPE 2: simple table (hm3 and hm4)
# Rows = models/preprocessing
# Columns = the 3 metrics (RMSE, RPD, R²)
# One independent scale per metric via ggnewscale
# ======================================================
make_heatmap_simple <- function(df_wide, row_var,
                                 title, subtitle, filename,
                                 width=10, height=7) {

  rmse_vals <- df_wide$RMSE_Test
  q33  <- quantile(rmse_vals,1/3,na.rm=TRUE)
  q66  <- quantile(rmse_vals,2/3,na.rm=TRUE)
  rmin <- min(rmse_vals,na.rm=TRUE)
  rmax <- max(rmse_vals,na.rm=TRUE)
  rpd_max <- max(df_wide$RPD_Test,na.rm=TRUE)

  lbl_rmse <- paste0("RMSE Test\n(", CONC_UNITS, ")")
  lbl_rpd  <- "RPD Test"
  lbl_r2   <- "R\u00b2 Test"

  df_long <- df_wide %>%
    rename(Row = !!sym(row_var)) %>%
    pivot_longer(cols=c(RMSE_Test, RPD_Test, R2_Test),
                 names_to="metric_id", values_to="Value") %>%
    mutate(Metric = factor(
      case_when(metric_id=="RMSE_Test" ~ lbl_rmse,
                metric_id=="RPD_Test"  ~ lbl_rpd,
                metric_id=="R2_Test"   ~ lbl_r2),
      levels=c(lbl_rmse, lbl_rpd, lbl_r2)
    ))

  d_rmse <- filter(df_long, metric_id=="RMSE_Test")
  d_rpd  <- filter(df_long, metric_id=="RPD_Test")
  d_r2   <- filter(df_long, metric_id=="R2_Test")

  p <- ggplot() +

    # ── RMSE ─────────────────────────────────────────────────────
    geom_tile(data=d_rmse, aes(x=Metric, y=Row, fill=Value),
              color="white", linewidth=0.7) +
    geom_text(data=d_rmse, aes(x=Metric, y=Row, label=round(Value,2)),
              color="black", fontface="bold", size=3.5) +
    scale_fill_gradientn(
      name    = paste0("RMSE\n(", CONC_UNITS, ")"),
      colours = c("darkgreen","yellowgreen","yellow","orange","red"),
      values  = rescale(c(rmin, q33, (q33+q66)/2, q66, rmax)),
      limits  = c(rmin, rmax),
      guide   = guide_colorbar(
        barheight=unit(3.5,"cm"), barwidth=unit(0.5,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=9,face="bold",color="black"),
        label.theme=element_text(size=8,color="black"), order=3)
    ) +
    new_scale_fill() +

    # ── RPD ──────────────────────────────────────────────────────
    geom_tile(data=d_rpd, aes(x=Metric, y=Row, fill=Value),
              color="white", linewidth=0.7) +
    geom_text(data=d_rpd, aes(x=Metric, y=Row, label=round(Value,2)),
              color="black", fontface="bold", size=3.5) +
    scale_fill_gradientn(
      name    = "RPD Test",
      colours = c("red","orange","yellow","yellowgreen","darkgreen"),
      values  = rescale(c(0, 1.0, 1.50, 2.50, rpd_max)),
      limits  = c(0, rpd_max),
      guide   = guide_colorbar(
        barheight=unit(3.5,"cm"), barwidth=unit(0.5,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=9,face="bold",color="black"),
        label.theme=element_text(size=8,color="black"), order=2)
    ) +
    new_scale_fill() +

    # ── R² ───────────────────────────────────────────────────────
    geom_tile(data=d_r2, aes(x=Metric, y=Row, fill=Value),
              color="white", linewidth=0.7) +
    geom_text(data=d_r2, aes(x=Metric, y=Row, label=round(Value,3)),
              color="black", fontface="bold", size=3.5) +
    scale_fill_gradientn(
      name    = "R\u00b2 Test",
      colours = c("red","orange","yellow","yellowgreen","darkgreen"),
      values  = rescale(c(0, 0.60, 0.80, 0.90, 1.0)),
      limits  = c(0, 1),
      guide   = guide_colorbar(
        barheight=unit(3.5,"cm"), barwidth=unit(0.5,"cm"),
        title.position="top", title.hjust=0.5,
        title.theme=element_text(size=9,face="bold",color="black"),
        label.theme=element_text(size=8,color="black"), order=1)
    ) +

    labs(title=title, subtitle=subtitle, x=NULL, y=NULL) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme_minimal(base_size=11) +
    theme(
      plot.title    = element_text(face="bold", hjust=0.5, size=13),
      plot.subtitle = element_text(hjust=0.5, size=9, color="gray40",
                                   margin=margin(b=6)),
      axis.text.x   = element_text(angle=0, hjust=0.5, size=10,
                                    face="bold", color="black"),
      axis.text.y   = element_text(size=10, color="black"),
      axis.ticks    = element_blank(),
      panel.grid    = element_blank(),
      legend.position   = "right",
      legend.box        = "vertical",
      legend.spacing.y  = unit(0.5,"cm"),
      legend.key.size   = unit(0.5,"cm"),
      plot.margin   = margin(8,5,8,5)
    )

  ggsave(filename, p, width=width, height=height, dpi=300,bg="white")
  cat("   OK Heatmap saved:", basename(filename), "\n")
  invisible(p)
}

# ======================================
# 7. Models
# ======================================
ctrl <- trainControl(method="cv", number=5, allowParallel=TRUE, verboseIter=FALSE)

# XGB is excluded from the models list and handled with a direct implementation
# (as in the classification pipeline) because caret + doParallel +
# xgboost raises "Error: Stopping" due to a fork conflict with OpenMP.
models <- list(
  PLS = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="pls",trControl=ctrl,tuneLength=15,metric="RMSE") },
  SVM = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="svmRadial",trControl=ctrl,
          preProcess=c("center","scale"),tuneLength=8,metric="RMSE") },
  Ridge = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="glmnet",trControl=ctrl,
          tuneGrid=expand.grid(alpha=0,lambda=10^seq(-2,2,length=40)),metric="RMSE") },
  Lasso = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="glmnet",trControl=ctrl,
          tuneGrid=expand.grid(alpha=1,lambda=10^seq(-4,0,length=40)),metric="RMSE") },
  ElasticNet = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="glmnet",trControl=ctrl,
          tuneGrid=expand.grid(alpha=seq(0.1,0.9,length=7),
                               lambda=10^seq(-3,1,length=15)),metric="RMSE") },
  RF = function(tr) { set.seed(123)
    train(Grupo~.,data=tr,method="rf",trControl=ctrl,
          tuneLength=6,ntree=300,metric="RMSE") }
)

# ── Direct XGB for REGRESSION (no caret, no doParallel) ──────────────
# Same strategy as in the classification pipeline:
# plain xgb.train() + manual 5-fold CV + nthread=1.
# Completely avoids the OpenMP/fork conflict that causes "Error: Stopping".
train_xgb_direct <- function(tr) {
  set.seed(123)

  y_num  <- as.numeric(tr$Grupo)
  X_mat  <- as.matrix(tr[, setdiff(names(tr), "Grupo")])
  dtrain <- xgb.DMatrix(data = X_mat, label = y_num)

  grid <- expand.grid(
    nrounds          = c(100, 200, 300),
    max_depth        = c(2, 3, 4),
    eta              = c(0.03, 0.05, 0.10),
    gamma            = c(0, 1),
    colsample_bytree = 0.7,
    min_child_weight = 10,
    subsample        = 0.7,
    stringsAsFactors = FALSE
  )

  best_rmse   <- Inf
  best_params <- NULL
  nfolds      <- 5
  folds       <- createFolds(y_num, k = nfolds, list = TRUE)

  for (i in seq_len(nrow(grid))) {
    params <- list(
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      max_depth        = grid$max_depth[i],
      eta              = grid$eta[i],
      gamma            = grid$gamma[i],
      colsample_bytree = grid$colsample_bytree[i],
      min_child_weight = grid$min_child_weight[i],
      subsample        = grid$subsample[i],
      nthread          = 1,   # no internal xgboost parallelism
      verbosity        = 0
    )

    cv_rmses <- numeric(nfolds)
    for (f in seq_len(nfolds)) {
      val_idx  <- folds[[f]]
      tr_idx   <- setdiff(seq_len(nrow(X_mat)), val_idx)
      d_tr     <- xgb.DMatrix(data  = X_mat[tr_idx,  , drop=FALSE],
                               label = y_num[tr_idx])
      d_val    <- xgb.DMatrix(data  = X_mat[val_idx, , drop=FALSE],
                               label = y_num[val_idx])
      m_tmp    <- xgb.train(params = params, data = d_tr,
                             nrounds = grid$nrounds[i], verbose = 0)
      pred_val <- predict(m_tmp, d_val)
      cv_rmses[f] <- sqrt(mean((pred_val - y_num[val_idx])^2, na.rm=TRUE))
    }

    mean_rmse <- mean(cv_rmses, na.rm=TRUE)
    if (mean_rmse < best_rmse) {
      best_rmse   <- mean_rmse
      best_params <- list(params=params, nrounds=grid$nrounds[i])
    }
  }

  final_model <- xgb.train(
    params  = best_params$params,
    data    = dtrain,
    nrounds = best_params$nrounds,
    verbose = 0
  )

  list(model=final_model, cv_rmse=best_rmse, best_params=best_params)
}

# Prediction for direct XGB (regression)
predict_xgb_reg <- function(xgb_obj, newdata_df) {
  X_mat <- as.matrix(newdata_df[, setdiff(names(newdata_df), "Grupo")])
  dmat  <- xgb.DMatrix(data = X_mat)
  as.numeric(predict(xgb_obj$model, dmat))
}

# ======================================
# 8. Experimental setup
# ======================================
scatter_opts <- c("none","SNV","MSC")
deriv_opts   <- c(0,1,2)
boruta_opts  <- c("none","boruta")

results_all     <- list()
error_log       <- list()
boruta_vars_log <- list()   # log of the variables selected by Boruta
start_time <- Sys.time()

cat("\n>>> EXPERIMENT CONFIGURATION:\n")
cat("   Preprocessing methods:",length(scatter_opts)*length(deriv_opts),"\n")
cat("   Models (caret):",length(models),"+ XGB (direct) = ",length(models)+1,"total\n")
cat("   Boruta options:",length(boruta_opts),"\n")
cat("   TOTAL COMBINATIONS:",
    length(scatter_opts)*length(deriv_opts)*length(boruta_opts)*(length(models)+1),"\n")
cat("   Estimated time: ~1-3 hours\n\n")

# ======================================
# 9. Main loop
# ======================================
total_runs  <- length(scatter_opts)*length(deriv_opts)*length(boruta_opts)*(length(models)+1)
current_run <- 0

for (sc in scatter_opts) {
  for (dv in deriv_opts) {
    pp_name<-make_pp_name(sc, dv)
    cat("\n>>> Preprocessing:",pp_name,"\n")
    pp<-apply_preprocessing(X_train,X_test,y_train,y_test,sc,dv)
    if(!pp$success){cat("   ERROR:",pp$error,"\n")
      error_log[[length(error_log)+1]]<-list(step=pp_name,error=pp$error);next}
    plot_spectra(pp$Xtr_plot,pp$wavelengths,paste("Preprocessed spectra (uncentered):",pp_name),NULL,pp_name)

    for (sel in boruta_opts) {
      train_base<-data.frame(Grupo=pp$ytr,pp$Xtr)
      test_base <-data.frame(Grupo=pp$yte,pp$Xte)

      if(sel=="boruta"){
        cat("   Running Boruta...\n")
        bor_res<-tryCatch({
          bor<-Boruta(Grupo~.,train_base,maxRuns=500,doTrace=0)
          n_tentative<-sum(bor$finalDecision=="Tentative")
          if(n_tentative>0){
            cat("   ",n_tentative,"variable(s) Tentative -> resolving with TentativeRoughFix...\n")
            bor<-TentativeRoughFix(bor)
          }
          vars<-names(bor$finalDecision[bor$finalDecision=="Confirmed"])
          if(length(vars)==0){cat("   WARNING: no vars confirmed. Using all.\n")
            list(success=TRUE,use_all=TRUE,vars=colnames(pp$Xtr),n_tentative=n_tentative)
          }else{cat("   Boruta selected",length(vars),"variables\n")
            list(success=TRUE,use_all=FALSE,vars=vars,n_tentative=n_tentative)}
        },error=function(e)list(success=FALSE,error=as.character(e)))
        if(!bor_res$success){
          error_log[[length(error_log)+1]]<-list(
            step=paste(pp_name,"Boruta"),error=bor_res$error);next}
        # Save the log of selected variables for the report
        boruta_vars_log[[paste(pp_name,"boruta",sep="__")]] <- list(
          Preprocessing = pp_name,
          N_selected    = length(bor_res$vars),
          N_tentative   = bor_res$n_tentative,
          Variables     = paste(bor_res$vars, collapse=", "),
          use_all       = if(exists("use_all",bor_res)) bor_res$use_all else FALSE
        )
        if(!bor_res$use_all){
          train_base<-train_base[,c("Grupo",bor_res$vars)]
          test_base <-test_base[, c("Grupo",bor_res$vars)]
          plot_spectra(pp$Xtr_plot,pp$wavelengths,paste("Boruta vars:",pp_name),bor_res$vars,pp_name)
        }
      }

      # ── caret models (PLS, SVM, Ridge, Lasso, ElasticNet, RF) ──────
      for (m in names(models)) {
        current_run <- current_run + 1
        cat("   [",current_run,"/",total_runs,"] Model:",m,"\n")
        res <- tryCatch({
          model <- models[[m]](train_base)
          ev    <- calc_metrics(model, train_base, test_base)
          if (!ev$success) return(list(success=FALSE, error=ev$error))
          plot_predictions(model, train_base, test_base, m, pp_name, sel)
          list(success=TRUE, metrics=ev)
        }, error=function(e) list(success=FALSE, error=as.character(e)))
        if (res$success) {
          ev <- res$metrics
          results_all[[length(results_all)+1]] <- data.frame(
            Preprocessing=pp_name, Boruta=sel, Model=m, N_vars=ncol(train_base)-1,
            R2_Train=ev$R2_Train, RMSE_Train=ev$RMSE_Train, RPD_Train=ev$RPD_Train,
            R2_Test =ev$R2_Test,  RMSE_Test =ev$RMSE_Test,  RPD_Test =ev$RPD_Test,
            Hyperparameters = format_bestTune(model$bestTune))
        } else {
          cat("   ERROR:",res$error,"\n")
          error_log[[length(error_log)+1]] <- list(step=paste(pp_name,sel,m), error=res$error)
        }
      }

      # ── Direct XGBoost (xgb.train, no caret, no doParallel) ──────
      # Run OUTSIDE the parallel cluster to avoid the OpenMP/fork
      # conflict that causes "Error: Stopping" with xgbTree in caret.
      current_run <- current_run + 1
      cat("   [",current_run,"/",total_runs,"] Model: XGB (direct xgb.train)\n")

      xgb_res <- tryCatch({

        # 1. Train
        xgb_obj <- train_xgb_direct(train_base)

        # 2. Predictions
        pt <- predict_xgb_reg(xgb_obj, train_base)
        pe <- predict_xgb_reg(xgb_obj, test_base)
        yt <- train_base$Grupo
        ye <- test_base$Grupo

        # 3. Metrics
        r2_tr  <- caret::R2(pt, yt);   rmse_tr <- caret::RMSE(pt, yt)
        rpd_tr <- sd(yt) / rmse_tr
        r2_te  <- caret::R2(pe, ye);   rmse_te <- caret::RMSE(pe, ye)
        rpd_te <- sd(ye) / rmse_te

        # 4. Scatter plot — same style as plot_predictions()
        pd <- data.frame(
          Real      = c(yt, ye),
          Predicted = c(pt, pe),
          Type      = c(rep("Training", length(yt)), rep("Test", length(ye)))
        )
        mn <- min(c(pd$Real, pd$Predicted), na.rm=TRUE)
        mx <- max(c(pd$Real, pd$Predicted), na.rm=TRUE)

        p_sc <- ggplot(pd, aes(x=Real, y=Predicted, color=Type)) +
          geom_point(size=3, alpha=0.7) +
          geom_abline(intercept=0, slope=1, linetype="dashed",
                      color="black", linewidth=1) +
          scale_color_manual(values=c("Training"="steelblue","Test"="firebrick")) +
          labs(title=paste0("XGB - ", pp_name),
               subtitle=paste0("Boruta: ", sel),
               x="Actual Concentration", y="Predicted Concentration",
               color="Dataset") +
          annotate("text", x=mn, y=mx,
                   label=paste0("TRAINING\nR2=",  round(r2_tr,3),
                                "\nRMSE=", round(rmse_tr,3),
                                "\nRPD=",  round(rpd_tr,2)),
                   hjust=0, vjust=1, size=3.5, color="steelblue", fontface="bold") +
          annotate("text", x=mx, y=mn,
                   label=paste0("TEST\nR2=",  round(r2_te,3),
                                "\nRMSE=", round(rmse_te,3),
                                "\nRPD=",  round(rpd_te,2)),
                   hjust=1, vjust=0, size=3.5, color="firebrick", fontface="bold") +
          theme_minimal(base_size=12) +
          theme(plot.title    = element_text(face="bold", hjust=0.5),
                plot.subtitle = element_text(hjust=0.5),
                legend.position="top") +
          coord_fixed()

        dir.create("Results/ScatterPlots/XGB", showWarnings=FALSE, recursive=TRUE)
        fn_sc <- paste0("Results/ScatterPlots/XGB/",
                        gsub(" ","_", pp_name), "_", sel, ".png")
        ggsave(fn_sc, p_sc, width=8, height=8, dpi=300,bg="white")

        list(success=TRUE,
             R2_Train=r2_tr, RMSE_Train=rmse_tr, RPD_Train=rpd_tr,
             R2_Test =r2_te, RMSE_Test =rmse_te, RPD_Test =rpd_te,
             Hyperparameters = format_xgb_params(xgb_obj$best_params))

      }, error=function(e) list(success=FALSE, error=as.character(e)))

      if (xgb_res$success) {
        results_all[[length(results_all)+1]] <- data.frame(
          Preprocessing=pp_name, Boruta=sel, Model="XGB", N_vars=ncol(train_base)-1,
          R2_Train=xgb_res$R2_Train, RMSE_Train=xgb_res$RMSE_Train, RPD_Train=xgb_res$RPD_Train,
          R2_Test =xgb_res$R2_Test,  RMSE_Test =xgb_res$RMSE_Test,  RPD_Test =xgb_res$RPD_Test,
          Hyperparameters = xgb_res$Hyperparameters)
        cat("   XGB OK — CV RMSE:", round(xgb_res$RMSE_Test, 4), "\n")
      } else {
        cat("   XGB ERROR:", xgb_res$error, "\n")
        error_log[[length(error_log)+1]] <- list(
          step  = paste(pp_name, sel, "XGB"),
          error = xgb_res$error)
      }

    }  # end boruta loop
  }    # end deriv loop
}      # end scatter loop

# ======================================
# 10. Export results and plots
# ======================================
if(length(results_all)>0){
  results_df<-bind_rows(results_all)
  write.xlsx(results_df,"Results/Excel/Complete_Summary.xlsx",overwrite=TRUE)
  cat("\nOK Results exported:",nrow(results_df),"rows\n")

  cat("\n=== TOP 10 (R2 Test) ===\n")
  print(results_df %>% arrange(desc(R2_Test)) %>% head(10) %>%
        select(Model,Preprocessing,Boruta,N_vars,R2_Test,RMSE_Test,RPD_Test))

  # ────────────────────────────────────────────────────────────────
  # Export the optimized hyperparameters of the top 10 models
  # (by R2_Test) to Excel, together with their metrics.
  # ────────────────────────────────────────────────────────────────
  top10_hp <- results_df %>%
    arrange(desc(R2_Test)) %>%
    head(10) %>%
    mutate(Rank = row_number(), .before = 1) %>%
    select(Rank, Model, Preprocessing, Boruta, N_vars, Hyperparameters,
           R2_Test, RMSE_Test, RPD_Test, R2_Train, RMSE_Train, RPD_Train)
  write.xlsx(top10_hp, "Results/Excel/Top10_Hyperparameters.xlsx", overwrite = TRUE)
  cat("   OK Saved: Top10_Hyperparameters.xlsx\n")

  # ────────────────────────────────────────────────────────────────
  # Export the variables selected by Boruta to Excel
  #   (a) All Preprocessing x Boruta combinations, one row
  #       per combination, with the full list of selected
  #       wavenumbers -> to see the overall effect of Boruta.
  #   (b) Only the combination of the WINNING MODEL (lowest RMSE_Test
  #       among the rows with Boruta=="boruta"), in long format
  #       (one wavenumber per row) so that it can later be cross-referenced
  #       with the functional groups of EDTA.
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
      arrange(RMSE_Test) %>%
      dplyr::slice(1)

    if (nrow(best_boruta_row) == 1) {
      best_key <- paste(best_boruta_row$Preprocessing, "boruta", sep = "__")
      best_log <- boruta_vars_log[[best_key]]

      if (!is.null(best_log)) {
        best_vars_wn <- suppressWarnings(parse_wn(
          trimws(strsplit(best_log$Variables, ",")[[1]])))
        best_vars_wn <- sort(best_vars_wn[!is.na(best_vars_wn)])

        best_vars_df <- data.frame(
          Winning_Model   = best_boruta_row$Model,
          Preprocessing   = best_boruta_row$Preprocessing,
          N_vars_selected = length(best_vars_wn),
          Wavenumber_cm1  = best_vars_wn,
          R2_Test         = round(best_boruta_row$R2_Test, 3),
          RMSE_Test       = round(best_boruta_row$RMSE_Test, 4),
          RPD_Test        = round(best_boruta_row$RPD_Test, 2))

        write.xlsx(best_vars_df,
                   "Results/Excel/Boruta_Variables_Best_Model.xlsx",
                   overwrite = TRUE)
        cat("   OK Saved: Boruta_Variables_Best_Model.xlsx  (Model:",
            best_boruta_row$Model, "| Preprocessing:", best_boruta_row$Preprocessing,
            "|", length(best_vars_wn), "wavenumbers )\n")
      } else {
        cat("   WARNING: no Boruta variable log found for the winning combination (",
            best_key, "). Skipping Boruta_Variables_Best_Model.xlsx\n")
      }
    } else {
      cat("   WARNING: no rows with Boruta=='boruta' in results_df. Skipping Boruta_Variables_Best_Model.xlsx\n")
    }
  }

  subtitle_hm<-paste0(
    "R\u00b2: red<0.80, yellow 0.80\u20130.90, green\u22650.90  |  ",
    "RPD: red<1.50, yellow 1.50\u20132.50, green\u22652.50  |  ",
    "RMSE (",CONC_UNITS,"): relative terciles (green=lowest)")

  cat("\n>>> Generating heatmaps...\n")

  # HM1: facet — Model × Preprocessing
  hm1<-results_df %>% group_by(Model,Preprocessing) %>%
    summarise(R2_Test=mean(R2_Test,na.rm=TRUE),RMSE_Test=mean(RMSE_Test,na.rm=TRUE),
              RPD_Test=mean(RPD_Test,na.rm=TRUE),.groups="drop")
  make_heatmap_facet(hm1,"Preprocessing","Model",
    "Performance Metrics (Test) \u2014 Model vs Preprocessing",
    subtitle_hm,"Results/Heatmaps/Metrics_Model_Preprocessing.png",width=14,height=11)

  # HM2: facet — Model × Config (Preprocessing + Boruta)
  hm2<-results_df %>% mutate(Config=paste(Preprocessing,Boruta,sep="\n")) %>%
    group_by(Model,Config) %>%
    summarise(R2_Test=mean(R2_Test,na.rm=TRUE),RMSE_Test=mean(RMSE_Test,na.rm=TRUE),
              RPD_Test=mean(RPD_Test,na.rm=TRUE),.groups="drop")
  make_heatmap_facet(hm2,"Config","Model",
    "Performance Metrics (Test) \u2014 Model vs Preprocessing + Boruta",
    subtitle_hm,"Results/Heatmaps/Metrics_Model_Boruta.png",width=22,height=11)

  # HM3: simple — average by Model (rows=models, cols=3 metrics)
  hm3<-results_df %>% group_by(Model) %>%
    summarise(R2_Test=mean(R2_Test,na.rm=TRUE),RMSE_Test=mean(RMSE_Test,na.rm=TRUE),
              RPD_Test=mean(RPD_Test,na.rm=TRUE),.groups="drop")
  make_heatmap_simple(hm3,"Model",
    "Performance Metrics (Test) \u2014 Average by Model",
    "Average across all preprocessing methods and Boruta options",
    "Results/Heatmaps/Metrics_by_Model.png",width=9,height=6)

  # HM4: simple — average by Preprocessing (rows=preprocessing, cols=3 metrics)
  hm4<-results_df %>% group_by(Preprocessing) %>%
    summarise(R2_Test=mean(R2_Test,na.rm=TRUE),RMSE_Test=mean(RMSE_Test,na.rm=TRUE),
              RPD_Test=mean(RPD_Test,na.rm=TRUE),.groups="drop")
  make_heatmap_simple(hm4,"Preprocessing",
    "Performance Metrics (Test) \u2014 Average by Preprocessing",
    "Average across all models and Boruta options",
    "Results/Heatmaps/Metrics_by_Preprocessing.png",width=9,height=7)

  # HM5: RMSE only
  hm5<-results_df %>% group_by(Model,Preprocessing) %>%
    summarise(RMSE_Test=mean(RMSE_Test,na.rm=TRUE),.groups="drop")
  q33s<-quantile(hm5$RMSE_Test,1/3); q66s<-quantile(hm5$RMSE_Test,2/3)
  p_hm5<-ggplot(hm5,aes(x=Preprocessing,y=Model,fill=RMSE_Test))+
    geom_tile(color="white",linewidth=0.7)+
    geom_text(aes(label=round(RMSE_Test,2)),color="black",fontface="bold",size=3.2)+
    scale_fill_gradientn(name=paste0("RMSE\n(",CONC_UNITS,")"),
                         colours=c("darkgreen","yellowgreen","yellow","orange","red"),
                         values=rescale(c(min(hm5$RMSE_Test),q33s,(q33s+q66s)/2,q66s,max(hm5$RMSE_Test))),
                         guide=guide_colorbar(barheight=8,barwidth=1))+
    labs(title=paste0("RMSE Test (",CONC_UNITS,") \u2014 Model vs Preprocessing"),
         subtitle="Green=lowest (best) | Yellow=middle | Red=highest (worst)",
         x="Preprocessing",y="Model")+
    theme_minimal(base_size=12)+
    theme(plot.title=element_text(face="bold",hjust=0.5,size=13),
          plot.subtitle=element_text(hjust=0.5,size=9,color="gray40"),
          axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
  ggsave("Results/Heatmaps/RMSE_Model_Preprocessing.png",p_hm5,width=12,height=8,dpi=300,bg="white")
  cat("   OK RMSE-only heatmap saved\n")

  # Model Comparison
  preproc_pal<-setNames(
    c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
      "#A65628","#F781BF","#999999","#00CED1")[seq_len(length(unique(results_df$Preprocessing)))],
    unique(results_df$Preprocessing))
  p_comp<-ggplot(results_df,aes(x=Model,y=R2_Test,color=Preprocessing))+
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=0.90,ymax=Inf,   fill="green", alpha=0.04)+
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=0.80,ymax=0.90,  fill="orange",alpha=0.06)+
    annotate("rect",xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0.80,  fill="red",   alpha=0.04)+
    geom_hline(yintercept=0.90,linetype="dashed",color="darkgreen", linewidth=0.5,alpha=0.7)+
    geom_hline(yintercept=0.80,linetype="dashed",color="darkorange",linewidth=0.5,alpha=0.7)+
    geom_jitter(aes(shape=Boruta),width=0.25,size=2.5,alpha=0.75)+
    stat_summary(aes(group=interaction(Model,Preprocessing)),
                 fun=median,geom="crossbar",width=0.5,linewidth=0.6,fatten=2,show.legend=FALSE)+
    scale_color_manual(values=preproc_pal)+
    scale_shape_manual(values=c("none"=16,"boruta"=17),
                       labels=c("none"="Without Boruta","boruta"="With Boruta"))+
    labs(title="R\u00b2 Test by Model and Preprocessing",
         subtitle="Points = individual combinations  |  Bar = median  |  Shapes = Boruta",
         x="Model",y="R\u00b2 Test",color="Preprocessing",shape="Variable selection")+
    theme_minimal(base_size=12)+
    theme(plot.title=element_text(face="bold",hjust=0.5,size=14),
          plot.subtitle=element_text(hjust=0.5,size=9,color="gray40"),
          axis.text.x=element_text(angle=45,hjust=1),
          legend.position="right",panel.grid.major.x=element_blank())
  ggsave("Results/Other/Model_Comparison.png",p_comp,width=13,height=7,dpi=300,bg="white")
  cat("   OK Model comparison saved\n")

  # Metrics PCA
  cat("\n>>> Generating PCA of metrics...\n")
  pca_data  <-results_df %>%
    select(R2_Train,RMSE_Train,RPD_Train,R2_Test,RMSE_Test,RPD_Test) %>% scale()
  pca_result<-prcomp(pca_data,center=FALSE,scale.=FALSE)
  pca_scores<-as.data.frame(pca_result$x[,1:2])
  colnames(pca_scores)<-c("PC1","PC2")
  pca_scores$Model        <-results_df$Model
  pca_scores$Preprocessing<-results_df$Preprocessing
  pca_scores$Boruta       <-results_df$Boruta
  ve<-summary(pca_result)$importance[2,1:2]*100
  p_pca<-ggplot(pca_scores,aes(x=PC1,y=PC2,color=Model,shape=Boruta))+
    geom_point(size=3,alpha=0.7)+
    scale_shape_manual(values=c("none"=16,"boruta"=15),
                       labels=c("none"="Without Boruta","boruta"="With Boruta"))+
    labs(title="PCA of Performance Metrics (Train + Test)",
         subtitle="R\u00b2, RMSE and RPD from both training and test sets",
         x=paste0("PC1 (",round(ve[1],1),"%)"),y=paste0("PC2 (",round(ve[2],1),"%)"),
         color="Model",shape="Variable Selection")+
    theme_minimal(base_size=12)+
    theme(plot.title=element_text(face="bold",hjust=0.5,size=14),
          plot.subtitle=element_text(hjust=0.5,size=10),
          legend.position="right",panel.border=element_rect(color="gray50",fill=NA))
  ggsave("Results/Other/PCA_Metrics.png",p_pca,width=14,height=8,dpi=300,bg="white")
  ld<-as.data.frame(pca_result$rotation[,1:2]);ld$Variable<-rownames(ld)
  write.xlsx(list("Scores"=pca_scores,"Loadings"=ld,
                  "Variance_Explained"=data.frame(
                    Component=paste0("PC",seq_along(pca_result$sdev)),
                    Var_Explained=summary(pca_result)$importance[2,]*100,
                    Cumulative   =summary(pca_result)$importance[3,]*100)),
             "Results/Excel/PCA_Analysis.xlsx",overwrite=TRUE)

  # Ranking
  cat("\n>>> Generating ranking...\n")
  ranking<-results_df %>% arrange(RMSE_Test) %>% mutate(Rank=row_number()) %>%
    select(Rank,Model,Preprocessing,Boruta,N_vars,
           R2_Test,RMSE_Test,RPD_Test,R2_Train,RMSE_Train,RPD_Train)
  ms<-results_df %>% group_by(Model) %>%
    summarise(N=n(),RMSE_Mean=mean(RMSE_Test),RMSE_SD=sd(RMSE_Test),
              RMSE_Min=min(RMSE_Test),R2_Mean=mean(R2_Test),
              R2_Max=max(R2_Test),RPD_Mean=mean(RPD_Test),.groups="drop") %>% arrange(RMSE_Mean)
  ps<-results_df %>% group_by(Preprocessing) %>%
    summarise(N=n(),RMSE_Mean=mean(RMSE_Test),R2_Mean=mean(R2_Test),
              RPD_Mean=mean(RPD_Test),.groups="drop") %>% arrange(RMSE_Mean)
  bs<-results_df %>% group_by(Boruta) %>%
    summarise(N=n(),RMSE_Mean=mean(RMSE_Test),R2_Mean=mean(R2_Test),
              RPD_Mean=mean(RPD_Test),.groups="drop")
  cr<-results_df %>% group_by(Model,Preprocessing) %>%
    summarise(RMSE_Test=mean(RMSE_Test),R2_Test=mean(R2_Test),
              RPD_Test=mean(RPD_Test),.groups="drop")

  # ======================================
  # Paired t-test: effect of Boruta
  #   To assess the effect of Boruta, a formal statistical test is used
  #   instead of only comparing means. A paired t-test is
  #   used (not an independent one) because, for each Model x
  #   Preprocessing combination, the "Boruta" and "None" versions come from
  #   the SAME wine samples evaluated with and without variable
  #   selection -> they are paired, not independent observations.
  #   A one-tailed test (alternative) is used because the hypothesis is
  #   directional: "Boruta improves performance", not just "changes it".
  # ======================================
  cat("\n>>> Paired t-test: Boruta effect...\n")
  boruta_lv <- sort(unique(results_df$Boruta))
  bor_lab   <- boruta_lv[grepl("bor", boruta_lv, ignore.case=TRUE)]
  none_lab  <- boruta_lv[!grepl("bor", boruta_lv, ignore.case=TRUE)]

  if(length(boruta_lv)==2 && length(bor_lab)==1 && length(none_lab)==1){

    wide_metric <- function(metric_col){
      results_df %>%
        select(Model,Preprocessing,Boruta,all_of(metric_col)) %>%
        tidyr::pivot_wider(names_from=Boruta, values_from=all_of(metric_col))
    }
    w_rmse <- wide_metric("RMSE_Test")
    w_r2   <- wide_metric("R2_Test")
    w_rpd  <- wide_metric("RPD_Test")

    # Pairs with NA in either arm are dropped (combinations
    # that failed in one of the two scenarios), and the number of
    # complete pairs available for the test is recorded.
    ok_rmse <- complete.cases(w_rmse[,c(bor_lab,none_lab)])
    ok_r2   <- complete.cases(w_r2[,c(bor_lab,none_lab)])
    ok_rpd  <- complete.cases(w_rpd[,c(bor_lab,none_lab)])

    pt_rmse <- t.test(w_rmse[[bor_lab]][ok_rmse], w_rmse[[none_lab]][ok_rmse],
                       paired=TRUE, alternative="less")     # H1: RMSE(Boruta) < RMSE(None)
    pt_r2   <- t.test(w_r2[[bor_lab]][ok_r2],     w_r2[[none_lab]][ok_r2],
                       paired=TRUE, alternative="greater")  # H1: R2(Boruta)   > R2(None)
    pt_rpd  <- t.test(w_rpd[[bor_lab]][ok_rpd],   w_rpd[[none_lab]][ok_rpd],
                       paired=TRUE, alternative="greater")  # H1: RPD(Boruta)  > RPD(None)

    boruta_ttest_df <- data.frame(
      Metric           = c("RMSE_Test","R2_Test","RPD_Test"),
      N_pairs          = c(sum(ok_rmse), sum(ok_r2), sum(ok_rpd)),
      Mean_Boruta      = c(mean(w_rmse[[bor_lab]][ok_rmse]), mean(w_r2[[bor_lab]][ok_r2]),  mean(w_rpd[[bor_lab]][ok_rpd])),
      Mean_None        = c(mean(w_rmse[[none_lab]][ok_rmse]),mean(w_r2[[none_lab]][ok_r2]), mean(w_rpd[[none_lab]][ok_rpd])),
      Mean_Diff        = c(unname(pt_rmse$estimate), unname(pt_r2$estimate), unname(pt_rpd$estimate)),
      t_statistic      = c(unname(pt_rmse$statistic), unname(pt_r2$statistic), unname(pt_rpd$statistic)),
      df               = c(unname(pt_rmse$parameter), unname(pt_r2$parameter), unname(pt_rpd$parameter)),
      p_value_one_tail = c(pt_rmse$p.value, pt_r2$p.value, pt_rpd$p.value),
      H1               = c("RMSE lower with Boruta","R2 higher with Boruta","RPD higher with Boruta"),
      Significant_0.05 = c(pt_rmse$p.value, pt_r2$p.value, pt_rpd$p.value) < 0.05,
      stringsAsFactors = FALSE)

    # two-sided p: the one-tailed test only evaluates 'Boruta improves'; if Boruta
    # significantly WORSENS performance (one-tailed p close to 1), this two-sided p shows it.
    boruta_ttest_df$p_value_two_sided <- 2 * pt(-abs(boruta_ttest_df$t_statistic), boruta_ttest_df$df)
    write.xlsx(boruta_ttest_df, "Results/Excel/Boruta_Paired_Ttest.xlsx", overwrite=TRUE)
    cat("   OK Paired t-test exported to Results/Excel/Boruta_Paired_Ttest.xlsx\n")
    print(boruta_ttest_df)
  } else {
    boruta_ttest_df <- NULL
    cat("   Warning: could not automatically identify the 2 Boruta levels ('boruta'/'none'); check results_df$Boruta manually.\n")
  }

  write.xlsx(list("Complete_Ranking"=ranking,"Top_20"=head(ranking,20),
                  "Summary_By_Model"=ms,"Summary_By_Preprocessing"=ps,
                  "Summary_By_Boruta"=bs,"Cross_Table"=cr),
             "Results/Excel/RMSE_Test_Ranking.xlsx",overwrite=TRUE)
  cat("OK Ranking exported\n")

  top10_plot<-ranking %>% head(10) %>%
    mutate(Combo=reorder(paste(Model,Preprocessing,Boruta,sep=" | "),-RMSE_Test)) %>%
    ggplot(aes(x=Combo,y=RMSE_Test,fill=Model))+
    geom_col()+geom_text(aes(label=round(RMSE_Test,3)),vjust=-0.5,size=3)+
    labs(title=paste0("Top 10 Best Combinations (Lowest RMSE Test, ",CONC_UNITS,")"),
         x="",y=paste0("RMSE Test (",CONC_UNITS,")"))+
    theme_minimal()+
    theme(plot.title=element_text(face="bold",hjust=0.5),
          axis.text.x=element_text(angle=45,hjust=1,size=8),legend.position="bottom")
  ggsave("Results/Other/Top10_Combinations.png",top10_plot,width=14,height=8,dpi=300,bg="white")
  cat("OK All plots generated\n")

}else{cat("\nError: No results generated\n")}

if(length(error_log)>0){
  write.xlsx(bind_rows(lapply(error_log,as.data.frame)),
             "Results/Excel/Error_Log.xlsx",overwrite=TRUE)
  cat("\nWarning: Errors:",length(error_log),"\n")
}

# ======================================
# 11. Final summary
# ======================================
cat("\n",rep("=",70),"\n",sep="")
cat("PIPELINE COMPLETED\n"); cat(rep("=",70),"\n",sep="")

end_time<-Sys.time()
duration<-difftime(end_time,start_time,units="mins")

if(length(results_all)>0){
  cat("\n EXECUTION SUMMARY:\n")
  cat("  - Start:", format(start_time,"%Y-%m-%d %H:%M:%S"),"\n")
  cat("  - End:",   format(end_time,   "%Y-%m-%d %H:%M:%S"),"\n")
  cat("  - Duration:",round(duration,2),"minutes\n")
  cat("  - Samples after outlier removal:",nrow(X),"\n")
  cat("  - Outliers removed:",n_outliers,"(T2 AND Q)\n")
  cat("  - Concentration units:",CONC_UNITS,"\n")
  cat("  - Successful combinations:",nrow(results_df),"\n")
  best<-results_df %>% arrange(RMSE_Test) %>% head(1)
  cat("\n BEST COMBINATION:\n")
  cat("  - Model:",        best$Model,"\n")
  cat("  - Preprocessing:",best$Preprocessing,"\n")
  cat("  - Boruta:",       best$Boruta,"\n")
  cat("  - RMSE Test:",    round(best$RMSE_Test,4),CONC_UNITS,"\n")
  cat("  - R2 Test:",      round(best$R2_Test,4),"\n")
  cat("  - RPD Test:",     round(best$RPD_Test,2),"\n")
}

cat("\n FILES GENERATED:\n")
cat("   Outliers/Raw_Spectra_Outliers.png\n")
cat("   Outliers/PCA_PC1_vs_PC2.png\n")
cat("   Outliers/PCA_PC2_vs_PC3.png\n")
cat("   Outliers/Influence_Plot_T2_vs_Q.png\n")
cat("   Heatmaps/Metrics_Model_Preprocessing.png  (facet: 3 rows/model)\n")
cat("   Heatmaps/Metrics_Model_Boruta.png          (facet: 3 rows/model)\n")
cat("   Heatmaps/Metrics_by_Model.png              (rows=models, cols=metrics)\n")
cat("   Heatmaps/Metrics_by_Preprocessing.png      (rows=prep., cols=metrics)\n")
cat("   Heatmaps/RMSE_Model_Preprocessing.png\n")
cat("   Other/Model_Comparison.png\n")
cat("   Other/Top10_Combinations.png\n")
cat("   Other/PCA_Metrics.png\n")
cat("   Excel/Complete_Summary.xlsx\n")
cat("   Excel/RMSE_Test_Ranking.xlsx\n")
cat("   Excel/PCA_Analysis.xlsx\n")

stopCluster(cl); registerDoSEQ()

write.table(
  data.frame(Start=format(start_time,"%Y-%m-%d %H:%M:%S"),
             End  =format(end_time,   "%Y-%m-%d %H:%M:%S"),
             Duration_min=round(duration,2),Samples_clean=nrow(X),
             Outliers_removed=n_outliers,Criterion="T2 AND Q",
             Conc_units=CONC_UNITS,
             Combinations=ifelse(length(results_all)>0,nrow(results_df),0),
             Errors=length(error_log)),
  "Results/Execution_Info.txt",row.names=FALSE,quote=FALSE)

# ============================================================
# PDF REPORT  /  REGRESSION PIPELINE
# Includes: full spectra + Boruta variables marked
#          with translucent red vertical lines
# Runs after the main pipeline and uses its objects
# Requires: install.packages(c("grid","gridExtra","png"))
# ============================================================

library(grid)
library(gridExtra)
library(png)



# ── PALETTE ───────────────────────────────────────────────────
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

# ── BASIC HELPERS ──────────────────────────────────────────
insert_png <- function(path) {
  if (!file.exists(path))
    return(textGrob(paste0("[Missing: ", basename(path), "]"),
                    gp = gpar(col = "red", fontsize = 9)))
  rasterGrob(readPNG(path), interpolate = TRUE,
             width = unit(1,"npc"), height = unit(1,"npc"))
}

draw_page_number <- function(n) {
  grid.lines(x = c(.05,.95), y = c(.028,.028),
             gp = gpar(col = "gray80", lwd = .5))
  grid.text(paste("Page", n), x = unit(.5,"npc"), y = unit(.016,"npc"),
            gp = gpar(fontsize = 9, col = "gray50"))
}

draw_section_bar <- function(label, y_top = .955) {
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
  grid.lines(x = c(.05,.95), y = c(y,y), gp = gpar(col = col, lwd = .6))

clean_cap <- function(fname) {
  b <- tools::file_path_sans_ext(basename(fname))
  trimws(gsub("  +", " ", gsub("_", " ", b)))
}

# ── SPECTRA PLOT WITH BORUTA LINES ───────────────────────────
# Builds a ggplot of the full spectrum (training set)
# and overlays translucent red vertical lines at the
# wavenumbers selected by Boruta.
make_spectra_boruta_plot <- function(X_matrix, wl_vec,
                                     sel_vars   = NULL,
                                     title_str  = "",
                                     y_col      = NULL) {
  
  # Convert matrix to long format
  Xm  <- as.matrix(X_matrix)
  nc  <- ncol(Xm); nr  <- nrow(Xm)
  wl  <- if (length(wl_vec) == nc) wl_vec else seq_len(nc)
  
  df <- data.frame(
    Wavenumber = rep(wl, each = nr),
    Absorbance = as.vector(Xm),
    Sample     = rep(seq_len(nr), times = nc)
  )
  df <- df[complete.cases(df), ]
  
  # Line color: if a continuous response is available, use it
  use_color <- !is.null(y_col) && length(y_col) == nr
  
  p <- ggplot(df, aes(x = Wavenumber, y = Absorbance, group = Sample))
  
  if (use_color) {
    df$y_val <- rep(y_col, times = nc)
    p <- ggplot(df, aes(x = Wavenumber, y = Absorbance,
                        group = Sample, color = y_val)) +
      scale_color_gradientn(
        name    = CONC_UNITS,
        colours = c("#1a237e","#1565c0","#26c6da","#a5d6a7","#ffee58","#ef6c00","#b71c1c"),
        guide   = guide_colorbar(barheight = unit(3,"cm"), barwidth = unit(.4,"cm"),
                                 title.position = "top", title.hjust = .5))
  }
  
  p <- p +
    geom_line(alpha = .30, linewidth = .35) +
    labs(title    = title_str,
         x        = expression("Wavenumber (cm"^{-1}*")"),
         y        = "Absorbance / Preprocessed signal") +
    theme_minimal(base_size = 12) +
    theme(plot.title       = element_text(face = "bold", hjust = .5, size = 13),
          panel.grid.minor = element_blank(),
          panel.border     = element_rect(color = "gray70", fill = NA),
          legend.position  = if (use_color) "right" else "none")
  
  # Translucent red vertical lines for Boruta variables
  if (!is.null(sel_vars) && length(sel_vars) > 0) {
    # Extract numeric wavenumbers
    sw <- suppressWarnings(parse_wn(sel_vars))
    sw <- sw[!is.na(sw) & sw %in% wl]
    if (length(sw) > 0) {
      vl_df <- data.frame(xintercept = sw)
      p <- p +
        geom_vline(data = vl_df,
                   aes(xintercept = xintercept),
                   color     = "red",
                   alpha     = .35,          # translucent
                   linewidth = .9,
                   inherit.aes = FALSE) +
        annotate("text",
                 x     = min(wl, na.rm = TRUE),
                 y     = max(df$Absorbance, na.rm = TRUE),
                 label = paste0(length(sw), " variables selected (Boruta)"),
                 hjust = 0, vjust = 1,
                 color = "red", size = 3.8, fontface = "bold")
    }
  }
  p
}

# Saves the plot as a temporary PNG and returns the grob
spectra_grob <- function(X_matrix, wl_vec, sel_vars = NULL,
                         title_str = "", y_col = NULL) {
  tmp <- tempfile(fileext = ".png")
  p   <- make_spectra_boruta_plot(X_matrix, wl_vec, sel_vars, title_str, y_col)
  tryCatch(
    ggplot2::ggsave(tmp, p, width = 11, height = 5, dpi = 150, bg = "white"),
    error = function(e) NULL
  )
  insert_png(tmp)
}

# ── CONCEPTUAL DIAGRAMS ───────────────────────────────────────
box_r <- function(x, y, w=.14, h=.055, fill=COL_ACCENT, label="", sz=8.5) {
  grid.roundrect(x=unit(x,"npc"),y=unit(y,"npc"),
                 width=unit(w,"npc"),height=unit(h,"npc"),r=unit(4,"pt"),
                 gp=gpar(fill=fill,col="white",lwd=1))
  grid.text(label,x=unit(x,"npc"),y=unit(y,"npc"),
            gp=gpar(fontsize=sz,fontface="bold",col="white"))
}
arrowh <- function(x0,x1,y,col=COL_ACCENT)
  grid.lines(x=c(x0,x1),y=c(y,y),
             arrow=arrow(length=unit(5,"pt"),type="closed"),
             gp=gpar(col=col,lwd=1.5,fill=col))

draw_diagram_pls <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("PLS — Latent Variable Decomposition",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  box_r(.14,.68,w=.13,h=.10,fill=COL_LIGHT,label="X matrix\n(spectra)")
  box_r(.14,.44,w=.13,h=.08,fill="#c0392b",label="y vector\n(conc.)")
  box_r(.42,.68,w=.16,h=.10,fill=COL_ACCENT,label="T scores\n(X space)")
  box_r(.42,.44,w=.16,h=.08,fill="#e67e22",label="U scores\n(y space)")
  box_r(.74,.56,w=.14,h=.08,fill="#27ae60",label="y hat\n(predicted)")
  arrowh(.21,.34,.68,COL_ACCENT); arrowh(.21,.34,.44,"#e67e22")
  arrowh(.50,.67,.68,COL_ACCENT); arrowh(.50,.67,.44,"#e67e22")
  grid.lines(x=c(.74,.74),y=c(.60,.68),gp=gpar(col=COL_ACCENT,lwd=1.2,lty="dashed"))
  grid.lines(x=c(.74,.74),y=c(.52,.44),gp=gpar(col="#e67e22",lwd=1.2,lty="dashed"))
  grid.text("ncomp: 1\u201315",x=.42,y=.28,gp=gpar(fontsize=8,col=COL_LIGHT,fontface="italic"))
  grid.text("Maximize cov(T, U)",x=.42,y=.18,gp=gpar(fontsize=7.5,col="gray50"))
  popViewport()
}

draw_diagram_svm <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("SVM-RBF — Feature Space & Epsilon Tube",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  set.seed(42)
  px<-runif(22,.08,.92); py<-runif(22,.25,.82)
  grid.circle(x=px,y=py,r=unit(4,"pt"),
              gp=gpar(fill=COL_ACCENT,col="white",lwd=.5,alpha=.75))
  grid.lines(x=c(.10,.90),y=c(.60,.52),gp=gpar(col=COL_MID,lwd=2))
  grid.lines(x=c(.10,.90),y=c(.68,.60),gp=gpar(col=COL_MID,lwd=1,lty="dashed"))
  grid.lines(x=c(.10,.90),y=c(.52,.44),gp=gpar(col=COL_MID,lwd=1,lty="dashed"))
  grid.text("\u03b5-tube",x=.92,y=.56,just=c("left","center"),
            gp=gpar(fontsize=8,col=COL_MID,fontface="bold"))
  grid.text("C & \u03c3 auto-tuned",x=.20,y=.20,
            gp=gpar(fontsize=7.5,col=COL_ACCENT))
  grid.text("RBF kernel",x=.75,y=.20,
            gp=gpar(fontsize=7.5,col="#c0392b"))
  popViewport()
}

draw_diagram_glmnet <- function(vp, model="Ridge") {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  ttl <- if(model=="Ridge")      "Ridge — L2: all coefs shrink, none = 0"
  else if(model=="Lasso") "Lasso — L1: sparse, many coefs = 0"
  else                    "Elastic Net — L1+L2: sparse + stable"
  grid.text(ttl,x=.5,y=.93,gp=gpar(fontsize=8.5,fontface="bold",col=COL_MID))
  set.seed(7)
  n_coef<-20; raw<-rnorm(n_coef,0,1)
  shrink<-if(model=="Ridge") raw*.45
  else if(model=="Lasso"){v<-raw*.55;v[abs(v)<.35]<-0;v}
  else {v<-raw*.50;v[abs(v)<.20]<-0;v*.85}
  bw<-.033; xs<-.08; yz<-.52; ys<-.18
  for(i in seq_along(shrink)){
    xc<-xs+(i-1)*(bw+.004)
    h<-abs(shrink[i])*ys
    yc<-if(shrink[i]>=0) yz+h/2 else yz-h/2
    clr<-if(shrink[i]==0) "gray80" else COL_ACCENT
    grid.rect(x=unit(xc,"npc"),y=unit(yc,"npc"),
              width=unit(bw,"npc"),height=unit(max(h,.003),"npc"),
              gp=gpar(fill=clr,col=NA))
  }
  grid.lines(x=c(.05,.95),y=c(yz,yz),gp=gpar(col="gray60",lwd=.8))
  grid.text("Coefficients",x=.50,y=yz-.22,gp=gpar(fontsize=7.5,col="gray50"))
  zp<-round(mean(shrink==0)*100)
  grid.text(paste0(zp,"% coefs = 0"),x=.70,y=.22,
            gp=gpar(fontsize=8,col=if(zp>0)"#c0392b" else COL_ACCENT,fontface="bold"))
  popViewport()
}

draw_diagram_rf <- function(vp) {
  pushViewport(vp)
  grid.rect(gp=gpar(fill="#f5f8fc",col="#dde6f0",lwd=.8))
  grid.text("Random Forest — Bootstrap Ensemble",
            x=.5,y=.93,gp=gpar(fontsize=9,fontface="bold",col=COL_MID))
  tx<-c(.18,.50,.82); fills<-c("#2980b9","#27ae60","#8e44ad")
  for(i in 1:3){
    grid.lines(x=c(tx[i],tx[i]),y=c(.28,.42),gp=gpar(col=fills[i],lwd=2))
    grid.lines(x=c(tx[i]-.10,tx[i],tx[i]+.10),y=c(.58,.42,.58),
               gp=gpar(col=fills[i],lwd=1.5))
    grid.lines(x=c(tx[i]-.06,tx[i]-.10),y=c(.68,.58),gp=gpar(col=fills[i],lwd=1.2))
    grid.lines(x=c(tx[i]+.06,tx[i]+.10),y=c(.68,.58),gp=gpar(col=fills[i],lwd=1.2))
    grid.circle(x=tx[i],y=.72,r=unit(12,"pt"),gp=gpar(fill=fills[i],col=NA,alpha=.85))
    grid.text(paste0("Tree ",i),x=tx[i],y=.20,
              gp=gpar(fontsize=7.5,col=fills[i],fontface="bold"))
    grid.text("Bootstrap",x=tx[i],y=.12,gp=gpar(fontsize=6.5,col="gray50"))
  }
  grid.lines(x=c(.18,.50,.82),y=c(.28,.22,.28),gp=gpar(col=COL_DARK,lwd=1.2,lty="dashed"))
  box_r(.50,.10,w=.18,h=.07,fill=COL_DARK,label="Average\nprediction",sz=7.5)
  popViewport()
}

draw_diagram_xgb <- function(vp) {
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
    grid.text(paste0("T",i),x=xs[i],y=.65,gp=gpar(fontsize=8,fontface="bold",col="white"))
    grid.text(paste0("r",i),x=xs[i],y=.44,gp=gpar(fontsize=6.5,col=fills[i]))
    if(i<4) arrowh(xs[i]+.05,xs[i+1]-.05,.65,COL_DARK)
  }
  grid.lines(x=c(.12,.82),y=c(.30,.30),gp=gpar(col=COL_DARK,lwd=1,lty="dashed"))
  box_r(.50,.18,w=.24,h=.07,fill=COL_DARK,label="F(x) = \u03a3 T_i",sz=7.5)
  grid.text("nthread=1 | xgb.train direct (no caret wrapper)",
            x=.50,y=.06,gp=gpar(fontsize=7,col="gray50",fontface="italic"))
  popViewport()
}

# ══════════════════════════════════════════════════════════════
# SAVE BORUTA SPECTRA OVER THE RAW SPECTRUM
# Only ONE PNG is generated per preprocessing: the original spectrum
# (post-outlier removal, without any preprocessing) with translucent red
# lines at the wavenumbers confirmed by Boruta.
# Generated BEFORE opening the PDF so that graphics devices are not mixed.
# ══════════════════════════════════════════════════════════════
cat(">>> Generating Boruta spectra overlay plots (raw spectra)...\n")
dir.create("Results/Spectra/Boruta_Overlays", showWarnings=FALSE, recursive=TRUE)

if (!exists("boruta_vars_log")) boruta_vars_log <- list()
y_col_train <- if(exists("y_train")) as.numeric(y_train) else NULL

# Complete raw spectrum — overview without marks (section cover)
raw_spectra_grob_path <- "Results/Spectra/Boruta_Overlays/_raw_all_samples.png"
if (exists("X_train") && exists("wavelengths")) {
  p_raw_all <- make_spectra_boruta_plot(
    X_matrix  = X_train,
    wl_vec    = wavelengths,
    sel_vars  = NULL,
    title_str = paste0("Raw Spectra — Training set (",
                       nrow(X_train)," samples, post-outlier removal)"),
    y_col     = y_col_train
  )
  tryCatch(ggsave(raw_spectra_grob_path, p_raw_all, width=12, height=5, dpi=150,bg="white"),
           error=function(e) NULL)
  cat("   OK Raw overview saved\n")
}

# One PNG per preprocessing: RAW spectrum + Boruta lines
boruta_png_map <- list()

for (sc in scatter_opts) {
  for (dv in deriv_opts) {
    pp_name_local <- make_pp_name(sc, dv)
    
    bor_key <- paste(pp_name_local, "boruta", sep="__")
    
    # Variables confirmed by Boruta for this preprocessing
    sel_vars_local <- NULL
    if (bor_key %in% names(boruta_vars_log)) {
      bvl <- boruta_vars_log[[bor_key]]
      if (nchar(bvl$Variables) > 0 && !isTRUE(bvl$use_all))
        sel_vars_local <- trimws(strsplit(bvl$Variables, ",")[[1]])
    }
    
    # PNG: raw spectrum (X_train) + red lines at the Boruta positions
    path_raw_bor <- paste0("Results/Spectra/Boruta_Overlays/",
                           gsub(" ","_", pp_name_local), "_raw_boruta.png")
    
    p_raw_bor <- make_spectra_boruta_plot(
      X_matrix  = X_train,
      wl_vec    = wavelengths,
      sel_vars  = sel_vars_local,
      title_str = paste0("Raw spectra + Boruta selection — ", pp_name_local),
      y_col     = y_col_train
    )
    tryCatch(ggsave(path_raw_bor, p_raw_bor, width=12, height=5, dpi=150,bg="white"),
             error=function(e) NULL)
    
    boruta_png_map[[pp_name_local]] <- list(
      raw_boruta = path_raw_bor,
      n_vars     = if(!is.null(sel_vars_local)) length(sel_vars_local) else NA
    )
    cat("   OK", pp_name_local,
        if(!is.null(sel_vars_local)) paste0("(", length(sel_vars_local), " vars)") else "(all vars)",
        "\n")
  }
}
cat(">>> Done.\n\n")

# ══════════════════════════════════════════════════════════════
# OPEN PDF
# ══════════════════════════════════════════════════════════════
pdf_path <- "Results/Scan_Regression_Report.pdf"
pdf(pdf_path, width=11, height=8.5, paper="USr")
pg <- 0L

# ══════════════════════════════════════════════════════════════
# P1 — COVER PAGE
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
grid.rect(x=0,y=1,width=1,height=.52,just=c("left","top"),
          gp=gpar(fill=COL_DARK,col=NA))
grid.text("FT-MIR Spectroscopy",x=.5,y=.90,
          gp=gpar(fontsize=26,fontface="bold",col="white"))
grid.text("Pipeline for Regression Screening",x=.5,y=.81,
          gp=gpar(fontsize=19,fontface="plain",col=COL_SILVER))
grid.text("Machine Learning-Based Quantitative Analysis Report",x=.5,y=.73,
          gp=gpar(fontsize=12,fontface="italic",col="#c8d8e8"))
grid.lines(x=c(.10,.90),y=c(.685,.685),gp=gpar(col=COL_SILVER,lwd=1.2))

grid.rect(x=unit(.06,"npc"),y=unit(.08,"npc"),
          width=unit(.52,"npc"),height=unit(.37,"npc"),
          just=c("left","bottom"),gp=gpar(fill="white",col=COL_MID,lwd=1.2))
grid.rect(x=unit(.06,"npc"),y=unit(.45,"npc"),
          width=unit(.52,"npc"),height=unit(.028,"npc"),
          just=c("left","top"),gp=gpar(fill=COL_MID,col=NA))
grid.text("Study Information",x=unit(.32,"npc"),y=unit(.436,"npc"),
          gp=gpar(fontsize=10,fontface="bold",col="white"))

REPORT_MATRIX     <- ""
REPORT_ANALYTE    <- sheet_name   # filled in automatically with the analyzed sheet
REPORT_INSTRUMENT <- ""
REPORT_N_SAMPLES  <- if(exists("X_raw")) as.character(nrow(X_raw)) else ""

fields <- list(
  list(label="Matrix:",                 val=REPORT_MATRIX),
  list(label="Analyte:",                val=REPORT_ANALYTE),
  list(label="Instrument:",             val=REPORT_INSTRUMENT),
  list(label="Number of total samples:",val=REPORT_N_SAMPLES)
)
for(i in seq_along(fields)){
  yy <- .405-(i-1)*.075
  grid.text(fields[[i]]$label,x=unit(.10,"npc"),y=unit(yy,"npc"),
            just=c("left","center"),gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
  grid.lines(x=c(.30,.555),y=c(yy-.012,yy-.012),
             gp=gpar(col="gray70",lwd=.7,lty="dotted"))
  if(nchar(fields[[i]]$val)>0)
    grid.text(fields[[i]]$val,x=unit(.31,"npc"),y=unit(yy,"npc"),
              just=c("left","center"),gp=gpar(fontsize=10,col="#333333"))
}
grid.text(paste("Generated:",format(Sys.time(),"%Y-%m-%d %H:%M")),
          x=unit(.94,"npc"),y=unit(.12,"npc"),
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
  list(n="1",  t="Outlier Detection & Removal",                         p=3),
  list(n="2",  t="Spectral Preprocessing Methods",                      p=5),
  list(n="3",  t="Data Split & Cross-Validation Strategy",              p=6),
  list(n="4",  t="Machine Learning Algorithms & Hyperparameters",       p=7),
  list(n="4",  t="  \u2514 Hyperparameter Summary Table (all models)",  p=14),
  list(n="5",  t="Boruta Variable Selection",                           p=15),
  list(n="6",  t="Boruta Variables on Raw Spectra",                     p=17),
  list(n="7",  t="Performance Results — Summary Tables",                p=28),
  list(n="8",  t="Heatmaps — Model vs Preprocessing",                   p=30),
  list(n="9",  t="Scatter Plots — by Algorithm",                        p=37)
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
  yy <- yy-.058
}
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P3 — OUTLIERS text
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("1. Outlier Detection & Removal")

y_cur <- draw_body(paste0(
  "Outlier detection was performed on the raw spectral matrix (",
  nrow(X_raw)," samples \u00d7 ",ncol(X_raw)," variables) using a two-criterion ",
  "PCA-based approach. Data were mean-centered; the number of PCs retained was the ",
  "minimum needed to explain \u226595% of cumulative spectral variance (",
  n_pcs," PCs, ",round(var_cum[n_pcs]*100,1),"% variance)."), y=.87)

y_cur <- draw_body(paste0(
  "Two criteria applied simultaneously: (1) Hotelling T\u00b2 compared against a ",
  "chi-squared threshold at 99% confidence (T\u00b2 limit = ",round(T2_lim,2),
  "). (2) Q residuals exceeding mean + 3 SD (Q limit = ",round(Q_lim,4),
  "). Only samples exceeding BOTH thresholds simultaneously (\"severe\" outliers) removed. ",
  n_outliers," sample(s) removed \u2192 ",sum(!flag_outlier)," retained."),
  y=y_cur-.025)

grid.rect(x=unit(.05,"npc"),y=unit(y_cur-.075,"npc"),
          width=unit(.90,"npc"),height=unit(.048,"npc"),
          just=c("left","bottom"),gp=gpar(fill="#eef4fb",col="#b0c8e0",lwd=.8))
grid.text(
  paste0("T\u00b2 flagged: ",sum(flag_T2),
         "   |   Q flagged: ",sum(flag_Q),
         "   |   Removed (T\u00b2 AND Q): ",n_outliers," of ",nrow(X_raw)),
  x=.5,y=y_cur-.052,
  gp=gpar(fontsize=9.5,fontface="bold",col=COL_MID))

grid.draw(editGrob(insert_png("Results/Outliers/Raw_Spectra_Outliers.png"),
                   vp=viewport(x=.5,y=.285,width=.90,height=.41)))
grid.text("Figure 1. Raw FTMIR spectra. Blue = retained; Red = outliers.",
          x=.5,y=.065,gp=gp_caption)
draw_page_number(pg)

# P4 PCA plots
grid.newpage(); pg <- pg + 1L
draw_section_bar("1. Outlier Detection & Removal (cont.)")
grid.draw(editGrob(insert_png("Results/Outliers/PCA_PC1_vs_PC2.png"),
                   vp=viewport(x=.25,y=.67,width=.46,height=.44)))
grid.draw(editGrob(insert_png("Results/Outliers/PCA_PC2_vs_PC3.png"),
                   vp=viewport(x=.75,y=.67,width=.46,height=.44)))
grid.text("Figure 2. PCA PC1 vs PC2.",x=.25,y=.44,gp=gp_caption)
grid.text("Figure 3. PCA PC2 vs PC3.",x=.75,y=.44,gp=gp_caption)
hrule(.425)
grid.draw(editGrob(insert_png("Results/Outliers/Influence_Plot_T2_vs_Q.png"),
                   vp=viewport(x=.5,y=.245,width=.55,height=.36)))
grid.text("Figure 4. Influence plot: T\u00b2 vs Q residuals. Dashed lines = thresholds.",
          x=.5,y=.055,gp=gp_caption)
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P5 — PREPROCESSING
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("2. Spectral Preprocessing Methods")

preproc_items <- list(
  list(name="Mean Centering  (applied to all combinations)",
       desc="Each spectral variable centered by subtracting the training-set column mean. Same mean applied to test set to prevent data leakage. Prerequisite for PLS and regularized methods."),
  list(name="Savitzky-Golay smoothing / derivatives  (0th, 1st, 2nd | poly=3, window=11)",
       desc=paste0("Computed on the complete, continuous spectrum (545 variables) before the Patz windows are cropped. ",
                   "0th = smoothed spectrum", if (SMOOTH_ALWAYS) " (applied to every combination, including those without derivative). " else " (NOT applied in this run: combinations without derivative are raw spectra). ",
                   "1st derivative removes additive baseline offsets. ",
                   "2nd derivative removes constant and linear baselines, enhances spectral resolution. ",
                   if (SCATTER_FIRST) "Order of operations: scatter correction (SNV/MSC) -> Savitzky-Golay -> Patz windows -> mean centering."
                   else "Order of operations: Savitzky-Golay -> scatter correction (SNV/MSC) -> Patz windows -> mean centering.")),
  list(name="Standard Normal Variate (SNV)",
       desc="Each spectrum scaled to zero mean and unit variance independently, correcting multiplicative scatter and path-length differences. Applied sample-wise; no reference spectrum required."),
  list(name="Multiplicative Scatter Correction (MSC)",
       desc="Linear regression of each spectrum against the training mean spectrum; additive and multiplicative scatter estimated and removed. Reference computed once from training set only.")
)
yy <- .865
for(pp_item in preproc_items){
  grid.text(paste0("\u25B6  ",pp_item$name),x=unit(.06,"npc"),y=unit(yy,"npc"),
            just=c("left","center"),gp=gp_sub)
  yy <- draw_body(pp_item$desc, y=yy-.032) - .022
}
hrule(yy-.01)
grid.text(paste0("Total combinations: ",length(scatter_opts)*length(deriv_opts),
                 "   (Scatter: none / SNV / MSC  \u00d7  Derivative: 0 / 1st / 2nd)"),
          x=.5,y=yy-.040,gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# P6 — DATA SPLIT & CV
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("3. Data Split & Cross-Validation Strategy")
n_tr<-if(exists("X_train")) nrow(X_train) else "N/A"
n_te<-if(exists("X_test"))  nrow(X_test)  else "N/A"

y_cur <- draw_body(paste0(
  "After outlier removal, the ",sum(!flag_outlier),"-sample dataset was partitioned into ",
  "training and test sets with a stratified 70/30 split. Stratification discretized y into ",
  "5 equally-spaced bins via cut(), then caret::createDataPartition(p=0.70) was applied ",
  "with set.seed(1234). This preserves the concentration distribution in both subsets."),
  y=.87)

y_cur <- draw_body(paste0(
  "All caret models used 5-fold CV (trainControl method='cv', number=5) on the training ",
  "set only for hyperparameter tuning. XGBoost was trained via xgb.train() with manual ",
  "5-fold CV and nthread=1. The test set was held out entirely for final evaluation. ",
  "Tuning metric: RMSE. Parallelization: ",num_cores," CPU cores (doParallel)."),
  y=y_cur-.025)

split_df<-data.frame(
  Parameter=c("Samples (post-outlier removal)","Training set","Test set",
              "Split ratio","Stratification","Random seed",
              "Cross-validation","Tuning metric"),
  Value=c(sum(!flag_outlier),n_tr,n_te,"70 / 30","Yes — 5 bins on y","1234",
          "5-fold CV (training only)","RMSE"),
  stringsAsFactors=FALSE)
tbl_s<-tableGrob(split_df,rows=NULL,
                 theme=ttheme_minimal(base_size=9.5,
                                      core   =list(fg_params=list(col=COL_DARK,hjust=0,x=.03)),
                                      colhead=list(fg_params=list(fontface="bold",col=COL_MID,hjust=0,x=.03))))
grid.draw(editGrob(tbl_s,
                   vp=viewport(x=.45,y=max(y_cur-.20,.28),width=.82,height=.38)))
draw_page_number(pg)

# ══════════════════════════════════════════════════════════════
# SECTION 4 — ALGORITHMS (one page each, with diagram)
# ══════════════════════════════════════════════════════════════
algo_list <- list(
  list(name="Partial Least Squares (PLS)", caret="method = 'pls'",
       desc=paste0("PLS is the reference method in chemometrics for high-dimensional, collinear spectral data. ",
                   "It decomposes both X (spectra) and y (concentrations) simultaneously into latent variables ",
                   "that maximize the covariance between score matrices T and U. Regression is performed in this ",
                   "compressed latent space. Particularly powerful when variables >> samples."),
       hyper=data.frame(Hyperparameter=c("ncomp (latent variables)"),
                        Range=c("1 to 15  (tuneLength = 15)"),stringsAsFactors=FALSE),
       diag=draw_diagram_pls),
  list(name="Support Vector Machine — Radial Basis (SVM-RBF)",
       caret="method = 'svmRadial', preProcess = c('center','scale')",
       desc=paste0("SVM-RBF maps inputs to a high-dimensional feature space via K(x,x') = exp(\u2212\u03c3||x\u2212x'||\u00b2) ",
                   "then finds a maximum-margin \u03b5-insensitive tube (epsilon-SVR). Samples outside the tube contribute ",
                   "to the loss. C controls the penalty for violations. \u03c3 estimated via sigest(). Robust to overfitting ",
                   "in high-dimensional spaces and handles non-linear spectral-property relationships."),
       hyper=data.frame(Hyperparameter=c("C (cost)","sigma (\u03c3, RBF width)"),
                        Range=c("Auto-estimated (tuneLength = 8)","Auto via sigest() heuristic"),
                        stringsAsFactors=FALSE),
       diag=draw_diagram_svm),
  list(name="Ridge Regression (L2 Regularization)",
       caret="method = 'glmnet', alpha = 0",
       desc=paste0("Ridge adds an L2 penalty (\u03bb\u03a3\u03b2j\u00b2) to the OLS objective, shrinking all coefficients ",
                   "toward zero proportionally. No coefficient is driven exactly to zero, so all spectral variables are ",
                   "retained. Particularly effective for collinear predictors: distributes weight across correlated bands."),
       hyper=data.frame(Hyperparameter=c("lambda (\u03bb, regularization)","alpha"),
                        Range=c("10^(-2) to 10^(2),  40 log-spaced values","Fixed = 0 (pure L2)"),
                        stringsAsFactors=FALSE),
       diag=function(vp) draw_diagram_glmnet(vp,"Ridge")),
  list(name="Lasso Regression (L1 Regularization)",
       caret="method = 'glmnet', alpha = 1",
       desc=paste0("Lasso adds an L1 penalty (\u03bb\u03a3|\u03b2j|) which induces exact sparsity: as \u03bb increases, ",
                   "coefficients are driven to exactly zero, performing automatic variable selection. Valuable for ",
                   "identifying the most informative spectral wavenumber regions and building parsimonious models."),
       hyper=data.frame(Hyperparameter=c("lambda (\u03bb, regularization)","alpha"),
                        Range=c("10^(-4) to 10^(0),  40 log-spaced values","Fixed = 1 (pure L1)"),
                        stringsAsFactors=FALSE),
       diag=function(vp) draw_diagram_glmnet(vp,"Lasso")),
  list(name="Elastic Net (L1 + L2 Regularization)",
       caret="method = 'glmnet', alpha in (0,1)",
       desc=paste0("Elastic Net combines L1 and L2 penalties: \u03bb[\u03b1\u03a3|\u03b2j|+(1-\u03b1)\u03a3\u03b2j\u00b2]. ",
                   "The mixing parameter \u03b1 interpolates between Lasso (\u03b1=1) and Ridge (\u03b1=0). ",
                   "Allows grouped variable selection: correlated spectral bands tend to be included or excluded together. ",
                   "7 \u03b1 levels \u00d7 15 \u03bb values = 105 combinations per preprocessing."),
       hyper=data.frame(Hyperparameter=c("alpha (\u03b1, L1/L2 mixing)","lambda (\u03bb, regularization)"),
                        Range=c("seq(0.1, 0.9)  [7 levels]","10^(-3) to 10^(1),  15 log-spaced values"),
                        stringsAsFactors=FALSE),
       diag=function(vp) draw_diagram_glmnet(vp,"ElasticNet")),
  list(name="Random Forest (RF)",
       caret="method = 'rf',  ntree = 300",
       desc=paste0("RF builds 300 CART trees, each on an independent bootstrap sample. At each node only mtry ",
                   "randomly chosen predictors are considered, decorrelating the trees. Predictions are averaged. ",
                   "RF provides implicit variable importance scores (% increase in MSE when a variable is permuted) ",
                   "and is robust to noise variables."),
       hyper=data.frame(Hyperparameter=c("mtry (variables per split)","ntree"),
                        Range=c("Auto: 6 candidate values  (tuneLength = 6)","Fixed = 300"),
                        stringsAsFactors=FALSE),
       diag=draw_diagram_rf),
  list(name="eXtreme Gradient Boosting (XGBoost)",
       caret="xgb.train() direct  |  objective = 'reg:squarederror'",
       desc=paste0("XGBoost builds trees sequentially, each fitting the negative gradient (pseudo-residuals) of the ",
                   "squared error loss. Contributions shrunk by eta. L1/L2 regularization, min_child_weight=10, and ",
                   "subsampling prevent overfitting. Implemented via xgb.train() directly (not caret) with nthread=1 ",
                   "to avoid the OpenMP/fork conflict causing 'Error: Stopping' inside doParallel workers."),
       hyper=data.frame(
         Hyperparameter=c("nrounds","max_depth","eta","gamma",
                          "colsample_bytree","min_child_weight","subsample"),
         Range=c("100, 200, 300","2, 3, 4","0.03, 0.05, 0.10","0, 1",
                 "Fixed = 0.70","Fixed = 10","Fixed = 0.70"),
         stringsAsFactors=FALSE),
       diag=draw_diagram_xgb)
)

for(al in algo_list){
  grid.newpage(); pg <- pg + 1L
  draw_section_bar(paste0("4. Algorithm: ",al$name))
  grid.text(paste0("caret / API: ",al$caret),
            x=unit(.06,"npc"),y=unit(.875,"npc"),
            just=c("left","center"),gp=gp_mono)
  hrule(.862)
  lines_d<-strwrap(al$desc,width=62)
  for(i in seq_along(lines_d))
    grid.text(lines_d[i],x=unit(.06,"npc"),y=unit(.835-(i-1)*.032,"npc"),
              just=c("left","top"),gp=gp_body)
  text_bot <- .835 - length(lines_d)*.032
  al$diag(viewport(x=.77,y=.67,width=.43,height=.38))
  hyp_y<-min(text_bot-.035,.44)
  grid.text("Hyperparameters evaluated:",
            x=unit(.06,"npc"),y=unit(hyp_y,"npc"),
            just=c("left","center"),gp=gp_sub)
  tbl_h<-tableGrob(al$hyper,rows=NULL,
                   theme=ttheme_minimal(base_size=9.5,
                                        core   =list(fg_params=list(col=COL_DARK,hjust=0,x=.02)),
                                        colhead=list(fg_params=list(fontface="bold",col=COL_MID,hjust=0,x=.02))))
  th<-nrow(al$hyper)*.052+.065
  grid.draw(editGrob(tbl_h,vp=viewport(x=.50,y=hyp_y-th/2-.03,width=.88,height=th)))
  draw_page_number(pg)
}

# ── Hyperparameter summary table (single page) ───────────
grid.newpage(); pg <- pg + 1L
draw_section_bar("4. Hyperparameter Summary — All Models")

draw_body(paste0(
  "The table below consolidates all hyperparameters evaluated across every model. ",
  "Fixed values were held constant for all combinations; tuned values were explored ",
  "via 5-fold cross-validation on the training set (metric: RMSE). ",
  "XGBoost used a full grid search with manual 5-fold CV via xgb.train()."),
  y=.875)

hyper_summary <- data.frame(
  Model = c(
    "PLS","PLS",
    "SVM-RBF","SVM-RBF",
    "Ridge","Ridge",
    "Lasso","Lasso",
    "ElasticNet","ElasticNet",
    "RF","RF",
    "XGBoost","XGBoost","XGBoost","XGBoost","XGBoost","XGBoost","XGBoost"
  ),
  Hyperparameter = c(
    "ncomp","[method]",
    "C (cost)","sigma (\u03c3)",
    "lambda (\u03bb)","alpha",
    "lambda (\u03bb)","alpha",
    "alpha (\u03b1)","lambda (\u03bb)",
    "mtry","ntree",
    "nrounds","max_depth","eta","gamma","colsample_bytree","min_child_weight","subsample"
  ),
  Values_Range = c(
    "1, 2, 3 \u2026 15  (15 values)","pls",
    "8 values auto-estimated from data","auto via sigest() heuristic",
    "10^\u22122 to 10\u00b2  (40 log-spaced values)","Fixed = 0  (pure L2)",
    "10^\u22124 to 10\u00b0  (40 log-spaced values)","Fixed = 1  (pure L1)",
    "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9  (7 values)","10^\u22123 to 10\u00b9  (15 log-spaced values)",
    "6 values auto-selected (tuneLength=6)","Fixed = 300",
    "100, 200, 300","2, 3, 4","0.03, 0.05, 0.10","0, 1","Fixed = 0.70","Fixed = 10","Fixed = 0.70"
  ),
  Type = c(
    "Tuned","Fixed",
    "Tuned","Tuned",
    "Tuned","Fixed",
    "Tuned","Fixed",
    "Tuned","Tuned",
    "Tuned","Fixed",
    "Tuned","Tuned","Tuned","Tuned","Fixed","Fixed","Fixed"
  ),
  stringsAsFactors = FALSE
)

# Colour rows by model
model_fills <- c(
  "PLS"="#dbeafe", "SVM-RBF"="#dcfce7", "Ridge"="#fef9c3",
  "Lasso"="#fce7f3", "ElasticNet"="#f3e8ff",
  "RF"="#ffedd5", "XGBoost"="#fee2e2"
)

tbl_hs <- tableGrob(hyper_summary, rows=NULL,
                    theme=ttheme_minimal(base_size=9,
                                         core    = list(fg_params=list(col=COL_DARK, hjust=0, x=.02),
                                                        bg_params=list(fill=c(
                                                          rep(model_fills["PLS"],       2),
                                                          rep(model_fills["SVM-RBF"],   2),
                                                          rep(model_fills["Ridge"],      2),
                                                          rep(model_fills["Lasso"],      2),
                                                          rep(model_fills["ElasticNet"], 2),
                                                          rep(model_fills["RF"],         2),
                                                          rep(model_fills["XGBoost"],    7)
                                                        ), col=NA)),
                                         colhead = list(fg_params=list(fontface="bold", col="white", hjust=.5, x=.5),
                                                        bg_params=list(fill=COL_MID, col=NA))))

grid.draw(editGrob(tbl_hs,
                   vp=viewport(x=.5, y=.45, width=.94, height=.62)))

# Color legend
leg_x <- .08
for(mn in names(model_fills)){
  grid.rect(x=unit(leg_x,"npc"), y=unit(.07,"npc"),
            width=unit(.09,"npc"), height=unit(.025,"npc"),
            gp=gpar(fill=model_fills[mn], col="gray60", lwd=.5))
  
  grid.text(mn,
            x=unit(leg_x + 0.045, "npc"),
            y=unit(.07,"npc"),
            just=c("left","center"),
            gp=gpar(fontsize=7.5, col=COL_DARK))
  
  leg_x <- leg_x + .13
}
draw_page_number(pg)
# ══════════════════════════════════════════════════════════════
grid.newpage(); pg <- pg + 1L
draw_section_bar("5. Boruta Variable Selection")

y_cur <- draw_body(paste0(
  "Boruta was run on the training set for each of the ",
  length(scatter_opts)*length(deriv_opts)," preprocessing combinations. ",
  "The algorithm wraps a Random Forest and iteratively confirms features whose importance ",
  "exceeds the maximum importance of randomly permuted shadow features (maxRuns=500). ",
  "Confirmed variables are retained; all others discarded before model training. ",
  "The table summarises the number confirmed and the % of the total spectral variables."),
  y=.87)

if(length(boruta_vars_log)>0){
  n_total_vars <- if(exists("X_raw")) ncol(X_raw) else NA
  bvl_df <- bind_rows(lapply(boruta_vars_log, function(x)
    data.frame(Preprocessing = x$Preprocessing,
               N_selected    = x$N_selected,
               N_total       = n_total_vars,
               Pct_selected  = paste0(round(x$N_selected/n_total_vars*100,1),"%"),
               stringsAsFactors=FALSE)))
  tbl_bv<-tableGrob(bvl_df,rows=NULL,
                    theme=ttheme_minimal(base_size=9,
                                         core   =list(fg_params=list(col=COL_DARK)),
                                         colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
  grid.draw(editGrob(tbl_bv,
                     vp=viewport(x=.5,y=max(y_cur-.18,.42),width=.75,height=.32)))
} else {
  grid.text("[boruta_vars_log not found — check pipeline loop]",
            x=.5,y=y_cur-.08,gp=gpar(fontsize=9,col="red",fontface="italic"))
}
draw_page_number(pg)

# ── Boruta variable list pages ────────────────────────────
if(length(boruta_vars_log)>0){
  chunks<-split(boruta_vars_log,ceiling(seq_along(boruta_vars_log)/7))
  for(ck in seq_along(chunks)){
    grid.newpage(); pg<-pg+1L
    draw_section_bar(paste0("5. Boruta — Variable Lists (page ",ck,")"))
    yy<-.88
    for(entry in chunks[[ck]]){
      grid.rect(x=unit(.05,"npc"),y=unit(yy,"npc"),
                width=unit(.90,"npc"),height=unit(.032,"npc"),
                just=c("left","top"),gp=gpar(fill="#eef4fb",col="#b0c8e0",lwd=.6))
      grid.text(paste0(entry$Preprocessing,
                       "  \u2014  ",entry$N_selected," variables confirmed"),
                x=unit(.07,"npc"),y=unit(yy-.016,"npc"),
                just=c("left","center"),
                gp=gpar(fontsize=9.5,fontface="bold",col=COL_MID))
      yy<-yy-.038
      var_text<-if(nchar(entry$Variables)>0) entry$Variables else "(all variables used)"
      vlines<-strwrap(var_text,width=115)
      for(vl in vlines){
        grid.text(vl,x=unit(.07,"npc"),y=unit(yy,"npc"),just=c("left","top"),
                  gp=gpar(fontsize=7.5,col="#444444",fontfamily="mono"))
        yy<-yy-.022
      }
      yy<-yy-.014
      if(yy<.10) break
    }
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# SECTION 6 — RAW SPECTRA WITH BORUTA VARIABLES MARKED
# One page per preprocessing: always on the original spectrum
# (post-outlier removal, without any mathematical transformation)
# ══════════════════════════════════════════════════════════════
if(length(boruta_png_map)>0){
  
  # Introduction page
  grid.newpage(); pg<-pg+1L
  draw_section_bar("6. Boruta Variable Selection — Highlighted on Raw Spectra")
  
  draw_body(paste0(
    "For each of the ",length(boruta_png_map)," preprocessing combinations, the figure shows ",
    "the raw training spectra (post-outlier removal, n=",n_tr," samples, no mathematical ",
    "transformation applied) with the variables confirmed by Boruta marked as translucent ",
    "red vertical lines (\u03b1=0.35, lwd=0.9)."),
    y=.86)
  
  draw_body(paste0(
    "Importantly, the Boruta selection was performed on the preprocessed version of the data ",
    "for each combination, but the selected wavenumber positions are always visualised here ",
    "on the original raw spectrum. This makes the chemical interpretation consistent: the ",
    "same physical absorption bands appear at the same x-axis positions regardless of ",
    "which preprocessing was used, enabling direct cross-comparison between combinations."),
    y=.74)
  
  # Interpretation box
  grid.rect(x=unit(.05,"npc"), y=unit(.26,"npc"),
            width=unit(.90,"npc"), height=unit(.27,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill="#f0f6ff", col=COL_ACCENT, lwd=1))
  grid.text("How to read these figures:",
            x=unit(.08,"npc"), y=unit(.51,"npc"),
            just=c("left","center"),
            gp=gpar(fontsize=10, fontface="bold", col=COL_MID))
  hints <- c(
    "\u25B6  Grey/coloured lines = raw training spectra (colour \u221d analyte concentration: blue\u2192red)",
    "\u25B6  Red vertical lines = wavenumber positions confirmed by Boruta for that preprocessing",
    "\u25B6  Line density reflects how many variables were selected (sparse = localised bands)",
    "\u25B6  Regions selected consistently across multiple preprocessings = chemically robust bands",
    "\u25B6  Each page title shows the preprocessing and the number of confirmed variables"
  )
  for(ii in seq_along(hints))
    grid.text(hints[ii],
              x=unit(.09,"npc"), y=unit(.475-(ii-1)*.040,"npc"),
              just=c("left","center"),
              gp=gpar(fontsize=9, col="#333333"))
  draw_page_number(pg)
  
  # Complete raw spectrum page (no marks — overview)
  grid.newpage(); pg<-pg+1L
  draw_section_bar("6. Raw Spectra Overview (no Boruta marks)")
  if(file.exists(raw_spectra_grob_path)){
    grid.draw(editGrob(insert_png(raw_spectra_grob_path),
                       vp=viewport(x=.5, y=.48, width=.93, height=.82)))
    grid.text(paste0("Figure. Raw training spectra — all ",n_tr," samples (post-outlier). ",
                     "Colour encodes analyte concentration (blue = low \u2192 red = high)."),
              x=.5, y=.062, gp=gp_caption)
  } else {
    draw_body("[Raw spectra overview PNG not found]", y=.5)
  }
  draw_page_number(pg)
  
  # One page per preprocessing: raw spectrum + Boruta lines
  fig_n <- 20L
  for(pp_nm in names(boruta_png_map)){
    grid.newpage(); pg<-pg+1L; fig_n<-fig_n+1L
    n_sel <- if(!is.na(boruta_png_map[[pp_nm]]$n_vars))
      boruta_png_map[[pp_nm]]$n_vars else "N/A"
    
    draw_section_bar(paste0("6. Boruta on Raw Spectra — ", pp_nm,
                            "  [", n_sel, " variables]"))
    
    grid.draw(editGrob(insert_png(boruta_png_map[[pp_nm]]$raw_boruta),
                       vp=viewport(x=.5, y=.48, width=.93, height=.82)))
    
    grid.text(
      paste0("Figure ",fig_n,". Raw training spectra (n=",n_tr,") with ",n_sel,
             " Boruta-confirmed wavenumbers marked in red (translucent). ",
             "Selection derived from preprocessed data: ",pp_nm,"."),
      x=.5, y=.062, gp=gp_caption)
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# RESULTS TABLES
# ══════════════════════════════════════════════════════════════
if(exists("results_df") && nrow(results_df)>0){
  
  grid.newpage(); pg<-pg+1L
  draw_section_bar("7. Performance Results — Top 20 Combinations")
  
  top20<-results_df %>% arrange(RMSE_Test) %>% head(20) %>%    mutate(Rank=row_number(),
                                                                      R2_Test   =round(R2_Test,3),   RMSE_Test=round(RMSE_Test,4),
                                                                      RPD_Test  =round(RPD_Test,2),  R2_Train =round(R2_Train,3),
                                                                      RMSE_Train=round(RMSE_Train,4)) %>%
    select(Rank,Model,Preprocessing,Boruta,N_vars,
           R2_Test,RMSE_Test,RPD_Test,R2_Train,RMSE_Train)
  
  tbl20<-tableGrob(top20,rows=NULL,
                   theme=ttheme_minimal(base_size=7.5,
                                        core   =list(fg_params=list(col=COL_DARK)),
                                        colhead=list(fg_params=list(fontface="bold",col=COL_MID,fontsize=8))))
  grid.draw(editGrob(tbl20,vp=viewport(x=.5,y=.47,width=.96,height=.76)))
  grid.text(paste0("Sorted by RMSE Test (ascending). Units: ",CONC_UNITS,"."),
            x=.5,y=.065,gp=gp_caption)
  draw_page_number(pg)
  
  grid.newpage(); pg<-pg+1L
  draw_section_bar("7. Performance Results — Summary by Model")
  
  ms_tbl<-results_df %>% group_by(Model) %>%
    summarise(N=n(),
              RMSE_Mean=round(mean(RMSE_Test,na.rm=TRUE),4),
              RMSE_SD  =round(sd(RMSE_Test,na.rm=TRUE),4),
              RMSE_Min =round(min(RMSE_Test,na.rm=TRUE),4),
              R2_Mean  =round(mean(R2_Test,na.rm=TRUE),3),
              R2_Max   =round(max(R2_Test,na.rm=TRUE),3),
              RPD_Mean =round(mean(RPD_Test,na.rm=TRUE),2),
              .groups="drop") %>% arrange(RMSE_Mean)
  
  tbl_ms<-tableGrob(ms_tbl,rows=NULL,
                    theme=ttheme_minimal(base_size=9.5,
                                         core   =list(fg_params=list(col=COL_DARK)),
                                         colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
  grid.draw(editGrob(tbl_ms,vp=viewport(x=.5,y=.77,width=.88,height=.28)))
  hrule(.61)
  grid.draw(editGrob(insert_png("Results/Other/Top10_Combinations.png"),
                     vp=viewport(x=.5,y=.355,width=.90,height=.44)))
  grid.text("Figure. Top 10 combinations ranked by RMSE Test (ascending).",
            x=.5,y=.065,gp=gp_caption)
  draw_page_number(pg)
}

# ══════════════════════════════════════════════════════════════
# HEATMAPS
# ══════════════════════════════════════════════════════════════
hm_list<-list(
  list(f="Results/Heatmaps/Metrics_Model_Preprocessing.png",
       c="Figure. R\u00b2, RMSE, RPD by Model vs Preprocessing (mean across Boruta options)."),
  list(f="Results/Heatmaps/Metrics_Model_Boruta.png",
       c="Figure. Metrics by Model vs Preprocessing + Boruta configuration."),
  list(f="Results/Heatmaps/Metrics_by_Model.png",
       c="Figure. Average performance metrics aggregated by model."),
  list(f="Results/Heatmaps/Metrics_by_Preprocessing.png",
       c="Figure. Average performance metrics aggregated by preprocessing method."),
  list(f="Results/Heatmaps/RMSE_Model_Preprocessing.png",
       c="Figure. RMSE Test heatmap: Model vs Preprocessing. Green = lowest (best); Red = highest."),
  list(f="Results/Other/Model_Comparison.png",
       c="Figure. R\u00b2 Test by model and preprocessing. Points = individual combinations; bar = median."),
  list(f="Results/Other/PCA_Metrics.png",
       c="Figure. PCA of performance metrics (R\u00b2, RMSE, RPD) from training and test sets.")
)
for(hm in hm_list){
  grid.newpage(); pg<-pg+1L
  draw_section_bar("8. Heatmaps — Model vs Preprocessing")
  grid.draw(editGrob(insert_png(hm$f),
                     vp=viewport(x=.5,y=.495,width=.92,height=.82)))
  grid.text(hm$c,x=.5,y=.062,gp=gp_caption)
  draw_page_number(pg)
}

# ══════════════════════════════════════════════════════════════
# SCATTER PLOTS (2×2 per page)
# ══════════════════════════════════════════════════════════════
for(m_name in c(names(models),"XGB")){
  scatter_dir<-paste0("Results/ScatterPlots/",m_name)
  if(!dir.exists(scatter_dir)) next
  png_files<-list.files(scatter_dir,pattern="\\.png$",full.names=TRUE)
  if(length(png_files)==0) next
  
  # Header page
  grid.newpage(); pg<-pg+1L
  draw_section_bar(paste0("9. Scatter Plots — ",m_name))
  draw_body(paste0(
    "Actual vs. Predicted concentration plots for all preprocessing \u00d7 Boruta combinations ",
    "evaluated with ",m_name,". Blue = Training set; Red = Test set. ",
    "Metrics shown in each panel: R\u00b2, RMSE (",CONC_UNITS,") and RPD."),
    y=.77)
  
  if(exists("results_df") && nrow(results_df)>0){
    best_m<-results_df %>% filter(Model==m_name) %>% arrange(RMSE_Test) %>% head(5) %>%
      mutate(R2_Test=round(R2_Test,3),RMSE_Test=round(RMSE_Test,4),RPD_Test=round(RPD_Test,2)) %>%
      select(Preprocessing,Boruta,N_vars,R2_Test,RMSE_Test,RPD_Test)
    tbl_b<-tableGrob(best_m,rows=NULL,
                     theme=ttheme_minimal(base_size=9,
                                          core   =list(fg_params=list(col=COL_DARK)),
                                          colhead=list(fg_params=list(fontface="bold",col=COL_MID))))
    grid.draw(editGrob(tbl_b,vp=viewport(x=.5,y=.47,width=.78,height=.28)))
    grid.text(paste0("Top 5 for ",m_name," sorted by RMSE Test."),
              x=.5,y=.305,gp=gp_caption)
  }
  draw_page_number(pg)
  
  # 2×2 grid pages
  n_plots<-length(png_files); ipp<-4L
  pos<-list(list(x=.25,y=.705),list(x=.75,y=.705),
            list(x=.25,y=.275),list(x=.75,y=.275))
  cap_y<-c(.488,.488,.055,.055)
  
  for(sp in seq_len(ceiling(n_plots/ipp))){
    grid.newpage(); pg<-pg+1L
    grid.text(paste0(m_name,"  \u2014  Page ",sp," of ",ceiling(n_plots/ipp)),
              x=.5,y=.974,gp=gpar(fontsize=10,fontface="bold",col=COL_MID))
    grid.lines(x=c(.04,.96),y=c(.955,.955),gp=gpar(col="gray75",lwd=.6))
    batch<-png_files[((sp-1)*ipp+1):min(sp*ipp,n_plots)]
    if(length(batch)>2) hrule(.495)
    for(bi in seq_along(batch)){
      grid.draw(editGrob(insert_png(batch[bi]),
                         vp=viewport(x=pos[[bi]]$x,y=pos[[bi]]$y,width=.48,height=.44)))
      grid.text(clean_cap(batch[bi]),x=pos[[bi]]$x,y=cap_y[bi],
                gp=gpar(fontsize=7.5,col="gray45",fontface="italic"))
    }
    draw_page_number(pg)
  }
}

# ══════════════════════════════════════════════════════════════
# CLOSE PDF
# ══════════════════════════════════════════════════════════════
dev.off()
cat("\n",rep("=",60),"\n",sep="")
cat("  REGRESSION PDF REPORT GENERATED\n")
cat("  Path:  ",pdf_path,"\n",sep="")
cat("  Pages: ",pg,"\n",sep="")
cat(rep("=",60),"\n\n",sep="")

# ══════════════════════════════════════════════════════════════
# Rename the results folder with the name of the analyzed
# sheet (e.g. "Results" -> "Results_POTASSIUM"), so that the
# results of another analyte are not overwritten in the next run.
# ══════════════════════════════════════════════════════════════
results_root <- paste0("Results_", sheet_name)
if (dir.exists(results_root)) {
  cat(">>> '",results_root,"' already existed (previous run for this same analyte) - overwriting.\n",sep="")
  unlink(results_root, recursive = TRUE)
}
if (dir.exists("Results")) {
  file.rename("Results", results_root)
  cat(">>> Results folder renamed to: ",results_root,"\n",sep="")
}

# ══════════════════════════════════════════════════════════════
# Append the analyte name to ALL generated files
# (e.g. "Boruta_Paired_Ttest.xlsx" -> "Boruta_Paired_Ttest_POTASSIUM.xlsx").
# This is automatic: the suffix comes from the name of the analyzed sheet (sheet_name).
# It is done at the end, once all files are closed, so as not to interfere
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
  cat(">>> Files renamed with suffix '_", tag, "': ", sum(ok), " of ", sum(todo), "\n", sep = "")
  if (sum(ok) < sum(todo))
    cat("    (files that could not be renamed are usually open in another program, e.g. Excel or a PDF viewer)\n")
}
