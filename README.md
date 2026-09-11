# Urban water demand modelling — Bellavista system, Quito (Ecuador)

Reproducible analysis code for the manuscript:

> *Nonlinear hydroclimatic and behavioural drivers of urban water demand:
> implications for peak-factor design in an Andean city.*
> Submitted to **Urban Water Journal**.

The pipeline models daily urban water demand for the Bellavista water supply
system using seven years of daily operational records (2007–2013), and
generates scenario-based projections to 2040.

---

## Quick start

```r
# 1. Open quito-water-demand.Rproj in RStudio, or set the working
#    directory to the folder containing run_forecasting.R
# 2. Install the dependencies (once)
install.packages(c("data.table", "forecast", "mgcv", "sandwich", "lmtest"))

# 3. Run
source("run_forecasting.R")
```

From a terminal:

```bash
Rscript run_forecasting.R
```

No paths need editing. The script resolves its own location and works on
Windows, macOS and Linux. Expect 20–30 minutes: the rolling-origin
cross-validation re-estimates ARIMA at 17 origins and the peak-factor
block runs 2,000 Gamma simulations.

---

## Repository layout

```
.
├── run_forecasting.R          Complete analysis pipeline
├── data/                      Input data (four CSV files)
├── outputs/
│   ├── figures/               Figures, internal naming
│   └── tables/                Main tables and their notes
├── results/                   Full result tables (CSV)
└── submission/                Generated: journal-ready package
    ├── figures/               Figures numbered as in the manuscript
    ├── tables/                Body tables
    ├── supplementary/         Supplementary tables S1–S6
    └── MANIFEST.txt           Contents and session information
```

`outputs/`, `results/` and `submission/` are produced by the script and are
not tracked by git.

---

## Data

| File | Variable | Source |
|---|---|---|
| `aguaQuito-dataset.csv` | Daily treated flow (L/s), calendar flags | EPMAPS |
| `C05-Bellavista_Precipitation-Daily.csv` | Daily rainfall (mm) | FONAG |
| `P09-Inaquito_INAMHI_Precipitation-Daily.csv` | Daily rainfall (mm) | INAMHI |
| `oni.csv` | Monthly Oceanic Niño Index | NOAA / NWS |

After quality control the analytical dataset contains **2,548 daily records**
(7 January 2007 – 31 December 2013). Three observations below the lower
boxplot fence, attributable to treatment-plant shutdowns rather than demand,
are excluded and reported by the script.

---

## What the pipeline does

1. **Reads and reconciles** the four sources, builds a composite
   precipitation series and imputes the few missing days by the median.
2. **Constructs predictors**: seven-day cumulative rainfall (P7) and its
   quadratic term, ONI interpolated to daily resolution, a dry-season
   indicator, structural trend and calendar controls.
3. **Estimates** univariate benchmarks (seasonal naive, ARIMA, ETS) and
   covariate specifications (M0 → MPOS) under GLM Gamma, OLS, GAM and
   SARIMAX, on a chronological 80/20 split.
4. **Evaluates** against two criteria declared in advance: predictive
   accuracy (rolling-origin cross-validation, MASE, horizons 1–365 days)
   and design fidelity (daily variability, peak-to-average ratio K).
5. **Refits** the preferred specification on the full observational record
   before generating any projection.
6. **Projects** to 2040 under rainfall and ENSO scenarios, with prediction
   intervals obtained by simulation from the fitted Gamma distribution.

---

## Figures and tables

| Manuscript | Internal name | Content |
|---|---|---|
| Figure 1 | — | Study-area map, prepared separately |
| Figure 2 | — | Methodological workflow, prepared separately |
| Figure 3 | `Figure4_ResidualDiagnostics` | Residual diagnostics, MPO |
| Figure 4 | `Figure2_ObservedVsPredicted` | Observed versus predicted |
| Figure 5 | `Figure6_HydrosocialDrivers` | Driver contributions |
| Figure 6 | `Figure3_RainfallDemand` | Response to antecedent rainfall |
| Figure 7 | `Figure5_ProjectedDemand` | Projections to 2040, two panels |
| Figure 8 | `Figure7_MASE_vs_Horizon` | Predictive accuracy, two panels |
| Figure 9 | `Figure_PeakFactor_K` | Annual peak factor K |
| Tables 1, 3, 5, 6 | — | Exported to `submission/tables/` |
| Tables S1–S6 | — | Exported to `submission/supplementary/` |

Internal names follow the order in which the figures are generated; the
export step renames them to the manuscript numbering in `submission/`.

Figures 1 and 2 must be placed in `submission/figures/` manually; the script
reports them as missing if absent.

---

## Reproducibility

Random seeds are fixed in the configuration block at the top of the script.
`MANIFEST.txt` records the R version and package versions of the run that
produced the package.

Tested with R 4.3 and R 4.5.

---

## Citation

Please cite the manuscript. A DOI will be added on acceptance.

## License

Code released under the MIT License (see `LICENSE`). The operational demand
data are provided by EPMAPS; please contact the authors regarding reuse.
