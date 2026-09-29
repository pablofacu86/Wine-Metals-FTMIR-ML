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
                              (POTASIO, MAGNESIO, CALCIO, HIERRO, COBRE)
pipelines/
  Pipeline_Scan_Regression.R
  Pipeline_Scan_Classification.R
reports/
  Scan_Regression_Report_Potasio.pdf      (pending)
  Scan_Regression_Report_Magnesio.pdf
  Scan_Regression_Report_Calcio.pdf
  Scan_Classification_Report_Hierro.pdf
  Scan_Classification_Report_Cobre.pdf
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
3. **Spectral preprocessing**: 9 combinations of Savitzky–Golay smoothing/
   derivatives (0th, 1st, 2nd) with Standard Normal Variate (SNV) or
   Multiplicative Scatter Correction (MSC), computed on the full spectrum,
   followed by extraction of the three spectral windows recommended by
   Patz et al. (965–1582, 1698–2006 and 2701–2971 cm⁻¹) and mean-centring.
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
   combinations whose test-set performance does not exceed the training
   performance by more than a fixed tolerance).
8. **Export of results**: Excel summaries, a multi-page PDF report, and
   diagnostic plots (spectra, PCA, heat maps, predicted-vs-measured plots,
   ROC curves, dendrograms, etc.), all written to a `Resultados_<ANALYTE>/`
   folder.

## How to run

1. Open `Pipeline_Scan_Regression.R` or `Pipeline_Scan_Classification.R` in
   R (≥ 4.3).
2. Edit the two lines near the top of the "Carga de datos" section:
   ```r
   data_path  <- "data/FINAL_DATA_SET.xlsx"   # path to the data file
   sheet_name <- "MAGNESIO"                   # analyte to run (see table below)
   ```
   | Pipeline | Valid `sheet_name` values |
   |---|---|
   | Scan Regression | `"POTASIO"`, `"MAGNESIO"`, `"CALCIO"` |
   | Scan Classification | `"HIERRO"`, `"COBRE"` |
3. Run the script. On completion, results are written to a folder named
   `Resultados_<sheet_name>/` (e.g. `Resultados_MAGNESIO/`), containing:
   - `Excel/` — `Complete_Summary.xlsx` (every combination evaluated),
     `Top10_Hyperparameters.xlsx`, `Boruta_Variables_All_Preprocessing.xlsx`,
     `Boruta_Variables_Best_Model.xlsx`, `Boruta_Paired_Ttest.xlsx`,
     `RMSE_Test_Ranking.xlsx` / ranking equivalent, `PCA_Analysis.xlsx`.
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

## Reports

The `reports/` folder contains the full PDF report generated by each
pipeline run (one per analyte), documenting the data set, preprocessing and
modelling choices, and the complete results for that analyte in detail.

## Citation

If you use this code or data, please cite:

> [Full citation of the article to be added once published.]

## License

[License to be defined by the authors — e.g. MIT for the code, CC-BY for the
data.]
