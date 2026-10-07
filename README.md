# FT-MIR metal estimation in wine — Scan Regression & Scan Classification

This repository contains the two R pipelines used to build and evaluate the
machine learning models described in *"Determination of Metals in Wines by
FT-MIR and machine learning tools"*, together with the raw FT-MIR spectra and
reference concentrations used to train them.

Both pipelines take the same kind of input (a wine FT-MIR spectrum, 902.57–
3000.84 cm⁻¹) and evaluate a systematic combination of spectral
preprocessing, variable selection and machine learning algorithms, in order
to compare their behaviour rather than to hand-pick a single "best" model.

- **Pipeline_Scan_Regression.R** — quantitative estimation of **potassium,
  magnesium and calcium** concentration (mg/L).
- **Pipeline_Scan_Classification.R** — binary classification of **iron and
  copper** relative to an enologically relevant threshold (10 mg/L and
  1 mg/L, respectively).

## Repository contents

```
data/
  FINAL_DATA_SET.xlsx        Raw spectra + reference concentrations, one sheet per analyte
                              (POTASSIUM, MAGNESIUM, CALCIUM, IRON, COPPER; the Spanish
                              names POTASIO, MAGNESIO, CALCIO, HIERRO, COBRE are also accepted)
pipelines/
  Pipeline_Scan_Regression.R
  Pipeline_Scan_Classification.R
results/
  Supplementary_Tables.xlsx  Results of every model evaluated, optimized hyperparameters,
                              Boruta-selected variables, paired t-tests and complete rankings
README.md
```

## Input data format

Each sheet of `FINAL_DATA_SET.xlsx` has the same layout:

| Column | Content |
|---|---|
| 1 | Sample identifier |
| 2 | Reference value: analyte concentration in mg/L (regression sheets) or class label (classification sheets) |
| 3+ | Absorbance at each of the 545 spectral variables, from 902.57 to 3000.84 cm⁻¹ in steps of ~3.86 cm⁻¹ |

The spectra are uploaded **without cropping**: each pipeline applies its own
spectral-window selection internally (see below), because smoothing and
derivatives must be computed on the continuous spectrum before any region is
removed.

## What each pipeline does

Both scripts follow the same overall workflow:

1. **Load data** from the Excel file and sheet configured at the top of the
   script (`data_path`, `sheet_name`).
2. **Outlier detection** on the raw spectra (restricted to the retained
   spectral windows) using a PCA-based Hotelling T² / Q-residual criterion;
   a sample is removed only if it exceeds **both** thresholds simultaneously.
3. **Spectral preprocessing**: 9 combinations (SG0, SG1, SG2, SNV + SG0,
   SNV + SG1, SNV + SG2, MSC + SG0, MSC + SG1, MSC + SG2). Scatter
   correction (SNV or MSC), when used, is applied first, followed by
   Savitzky–Golay (polynomial order 3, 11-point window) as smoothing (SG0)
   or as first (SG1) or second (SG2) derivative. Everything is computed on
   the full 545-variable spectrum, followed by extraction of the three
   spectral windows recommended by Patz et al. (965–1582, 1698–2006 and
   2701–2971 cm⁻¹) and mean-centring. The order and the smoothing can be
   changed with the switches `SCATTER_FIRST` and `SMOOTH_ALWAYS` at the top
   of each script.
4. **Variable selection** with the Boruta algorithm (500 iterations),
   evaluated both with and without selection.
5. **Train/test split**: stratified 70/30 split (by concentration bin for
   regression, by class for classification).
6. **Model training and 5-fold cross-validation** hyperparameter tuning for
   a panel of algorithms:
   - Regression: PLS, Ridge, Lasso, Elastic Net, SVM-RBF, Random Forest,
     XGBoost.
   - Classification: PLS-DA, LDA, SVM-RBF, Random Forest, CART, kNN, Naive
     Bayes, XGBoost.
7. **Evaluation** of every preprocessing × variable-selection × algorithm
   combination on the held-out test set, and identification of the best
   model per analyte (lowest test RMSE, or highest test AUC, among the
   combinations whose test-set R² or AUC does not exceed the training value
   by more than 0.05). In classification, the positive class is the one
   above the limit ("Higher than ..."), so sensitivity, specificity and F1
   refer to detecting samples over the limit.
8. **Export of results**: Excel summaries, a multi-page PDF report, and
   diagnostic plots (spectra, PCA, heat maps, predicted-vs-measured plots,
   ROC curves, dendrograms, etc.), all written to a `Results_<ANALYTE>/`
   folder.

## How to run

1. Open `Pipeline_Scan_Regression.R` or `Pipeline_Scan_Classification.R` in
   R (≥ 4.3).
2. Edit the two lines near the top of the data-loading section:
   ```r
   data_path  <- "data/FINAL_DATA_SET.xlsx"   # path to the data file
   sheet_name <- "MAGNESIUM"                  # analyte to run (see table below)
   ```
   | Pipeline | Valid `sheet_name` values |
   |---|---|
   | Scan Regression | `"POTASSIUM"`, `"MAGNESIUM"`, `"CALCIUM"` |
   | Scan Classification | `"IRON"`, `"COPPER"` |
3. Run the script. On completion, results are written to a folder named
   `Results_<sheet_name>/` (e.g. `Results_MAGNESIUM/`), containing:
   - `Excel/` — `Complete_Summary.xlsx` (every combination evaluated),
     `Top10_Hyperparameters.xlsx`, `Boruta_Variables_All_Preprocessing.xlsx`,
     `Boruta_Variables_Best_Model.xlsx`, `Boruta_Paired_Ttest.xlsx`,
     `RMSE_Test_Ranking.xlsx` / `Ranking.xlsx`, `PCA_Analysis.xlsx`. Every
     file name ends with the analyte name (e.g.
     `Complete_Summary_MAGNESIUM.xlsx`).
   - `Outliers/`, `Spectra/`, `Heatmaps/`, `ScatterPlots/` or
     `ConfusionMatrices/` / `ROC/`, `Other/` — diagnostic and results plots.
   - `Scan_Regression_Report.pdf` / `Scan_Classification_Report.pdf` — a
     single PDF summarising the whole run (data set, preprocessing,
     modelling choices, results and rankings).

### Requirements

R ≥ 4.3, with the following packages installed: `doParallel`, `dplyr`,
`tidyr`, `caret`, `pracma`, `pls`, `e1071`, `randomForest`, `openxlsx`,
`prospectr`, `Boruta`, `xgboost`, `ggplot2`, `glmnet`, `ggrepel`,
`ggnewscale`, `scales`, `pROC`, `MASS`, `rpart`, `naivebayes`, `Matrix`,
`grid`, `gridExtra`, `png`, `dendextend`, `circlize` (the last two are
installed automatically by the classification pipeline if missing, since
they are only needed for the circular dendrogram plot).

## Main results (version 1.1.0)

Best model of each analyte, selected as described in the article (lowest test
RMSE or highest test AUC among the combinations whose test R² or AUC does not
exceed the training value by more than 0.05):

| Analyte | Best model | Preprocessing | Boruta (variables) | Test-set performance |
|---|---|---|---|---|
| Potassium | Lasso | SG0 | no (310) | R² = 0.946; RMSE = 61.7 mg/L; RPD = 4.34 |
| Magnesium | PLS | SNV + SG1 | no (310) | R² = 0.836; RMSE = 7.02 mg/L; RPD = 2.42 |
| Calcium | XGBoost | SG1 | yes (59) | R² = 0.892; RMSE = 27.8 mg/L; RPD = 2.96 |
| Iron (10 mg/L) | SVM-RBF | MSC + SG0 | yes (46) | AUC = 0.985; accuracy = 0.915; κ = 0.816 |
| Copper (1 mg/L) | XGBoost | SNV + SG1 | yes (53) | AUC = 0.998; accuracy = 0.957; κ = 0.908 |

Several combinations per analyte perform almost as well as the selected one;
the 50 best models of each analyte are listed in the Supplementary Material of
the article and the complete rankings are in `results/Supplementary_Tables.xlsx`.
The results come from a single training/test partition (`set.seed(1234)`).

## Versions

The version described here (1.1.0, DOI [10.5281/zenodo.23199425](https://doi.org/10.5281/zenodo.23199425))
supersedes 1.0.1 (DOI 10.5281/zenodo.23033275).
Changes with respect to 1.0.1: scatter correction (SNV/MSC) is now applied before the
Savitzky–Golay filter, all preprocessing conditions include smoothing (SG0 = smoothing
only), the positive class in classification is the class above the limit, and all
results were regenerated. Version 1.0.1 used an earlier preprocessing order and should
not be used to reproduce the article.
The R version used to obtain the results was 4.3.1 (caret 6.0-94). Running
`sessionInfo()` at the end of a pipeline run records the exact package versions
of the local installation.

## Web app

The winning model of each analyte is available as a web application
(WineMetalScan): <https://github.com/pablofacu86/WineMetalScan>. The deployed models
were rebuilt from `data/FINAL_DATA_SET.xlsx` with the same
preprocessing and hyperparameters reported in the article.

## Citation

If you use this code or data, please cite the archived version of this repository:
https://doi.org/10.5281/zenodo.23199425
