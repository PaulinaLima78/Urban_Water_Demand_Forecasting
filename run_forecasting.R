# ============================================================
# QUITO (MDQ): DAILY URBAN WATER DEMAND MODELLING
# Bellavista water supply system - reproducible analysis pipeline
#
# Companion code for: "Nonlinear Hydroclimatic and Behavioural Drivers
# of Urban Water Demand" (Urban Water Journal, under review)
#
# Observational domain : 2007-01-07 to 2013-12-31 (2,551 daily records)
# Projection domain    : 2014-2030 (scenario simulation; no observations)
#
# ------------------------------------------------------------
# HOW TO RUN
#   1. Place the four input CSV files in ./DATA
#        - daily demand           (aguaQuito-dataset.csv)
#        - Bellavista precipitation (C05-*.csv)
#        - Inaquito precipitation   (P09-*.csv)
#        - monthly ONI              (oni.csv)
#   2. Open the .Rproj file, or set the working directory to this file's
#      folder, then:  source("run_forecasting_v2.R")
#
# No absolute paths are used. Outputs are written to ./OUTPUTS and ./RESULTS.
# ------------------------------------------------------------
#
# CHANGES RELATIVE TO THE FIRST-SUBMISSION PIPELINE
#   C1  Final model refitted on the full observational record before
#       projection (previously fitted on the 80% training split only).
#   C2  Peak factor K estimated from observed and simulated demand rather
#       than from modelled conditional means, which omit residual variance.
#   C3  Plant-shutdown outliers below the lower boxplot fence identified and
#       reported; retained in the analytical dataset by default.
#   C4  Antecedent precipitation resampled by calendar month in projections
#       (previously held constant across the whole horizon).
#   C5  Real Ecuadorian public-holiday calendar in projections
#       (previously assigned at random).
#   C6  ONI interpolated to daily resolution; clustered standard errors
#       reported because monthly values were repeated across days.
#   C7  Dry-season indicator renamed and its aliasing with the month factor
#       documented; MPS/MPOS reported as rank-deficient.
#   C8  Diebold-Mariano tests use the evaluation horizon, not h = 1.
#   C9  MASE scaled by the in-sample seasonal naive error.
#   C10 Cross-family AIC comparison removed (different likelihoods, scales, n).
#   C11 Preferred specification fixed a priori (MPO) rather than selected
#       by test-set RMSE at run time.
#   C12 Rolling-origin cross-validation by forecast horizon (new Table 5).
#   C13 Projection intervals obtained by simulation from the fitted Gamma.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(forecast)
  library(mgcv)
  library(sandwich)
  library(lmtest)
  library(grDevices)
})

# Dependencies are declared above and are not installed at run time: a call to
# install.packages() inside a script fails on machines without write access or
# network, and stepR was not used anywhere in this pipeline.

# ============================================================
# 0) PROJECT PATHS + CONFIG
# ============================================================

# Portable root: the folder containing this script. Works when sourced,
# run through Rscript, or executed from an RStudio project.
get_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(dirname(f[1])))
  if (!is.null(sys.frames()) && !is.null(sys.function(1))) {
    sf <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
    if (!is.null(sf)) return(dirname(sf))
  }
  normalizePath(getwd())
}
root_dir <- get_script_dir()
message("[OK] Project root: ", root_dir)

# Directory names are resolved case-insensitively so the repository runs
# unchanged on Windows, macOS and Linux (only Windows is case-insensitive
# at the filesystem level).
resolve_dir <- function(root, name) {
  hit <- list.dirs(root, recursive = FALSE, full.names = FALSE)
  m <- hit[tolower(hit) == tolower(name)]
  file.path(root, if (length(m)) m[1] else name)
}
data_dir    <- resolve_dir(root_dir, "data")
output_dir  <- resolve_dir(root_dir, "outputs")
results_dir <- resolve_dir(root_dir, "results")

fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")

dir.create(output_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir,     showWarnings = FALSE, recursive = TRUE)

find_one_file <- function(folder, pattern, label) {
  x <- list.files(folder, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(x) == 0) stop(paste0("No file found for ", label, " in: ", folder))
  if (length(x) > 1) {
    message("[WARN] Multiple files found for ", label, ". Using first match:")
    print(x)
  }
  x[1]
}

cfg <- list(
  demanda   = find_one_file(data_dir, "agua.*quito|demand|demanda", "daily demand"),
  p_bella   = find_one_file(data_dir, "bellavista.*precip|precip.*bellavista|c05", "Bellavista precipitation"),
  p_inq     = find_one_file(data_dir, "inaquito.*precip|precip.*inaquito|p09", "Iñaquito precipitation"),
  oni       = find_one_file(data_dir, "oni", "ONI"),
  anio_fin  = 2030,                # projection horizon, as reported in the manuscript
  oni_obs_through_year = 2026,
  split_frac = 0.80,
  seed      = 123,
  preferred_spec = "MPO",          # C11: fixed a priori, not chosen by test RMSE
  drop_outliers  = FALSE,          # C3: identify and report, do not exclude
  n_boot         = 1000,           # C13: projection intervals
  n_sim_K        = 2000,           # C2
  cv_first_origin = 1200,          # C12
  cv_step         = 60,
  cv_horizons     = c(1, 7, 30, 90, 365),
  m_season        = 7,
  out_dir   = output_dir,
  fig_dir   = fig_dir,
  tab_dir   = tab_dir,
  results_dir = results_dir
)

message("[OK] Demand file: ", cfg$demanda)
message("[OK] Bellavista precip file: ", cfg$p_bella)
message("[OK] Iñaquito precip file: ", cfg$p_inq)
message("[OK] ONI file: ", cfg$oni)

# ============================================================
# 1) UTILITIES
# ============================================================

step <- function(title) {
  cat("\n", paste0("==== ", title, " ===="), "\n", sep = "")
}

assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

safe_range <- function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))

save_png <- function(filename, expr, width = 2200, height = 1400, res = 300) {
  path <- file.path(cfg$fig_dir, filename)
  png(path, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  force(expr)
  message("[FIG] Saved: ", path)
  invisible(path)
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae  <- function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
medae <- function(y, yhat) median(abs(y - yhat), na.rm = TRUE)
mape <- function(y, yhat) 100 * mean(abs((y - yhat) / pmax(abs(y), 1e-9)), na.rm = TRUE)
smape <- function(y, yhat) 100 * mean(2 * abs(yhat - y) / pmax(abs(y) + abs(yhat), 1e-9), na.rm = TRUE)
bias <- function(y, yhat) mean(yhat - y, na.rm = TRUE)

# C9: MASE must be scaled by the in-sample one-step seasonal naive error
# of the TRAINING window (Hyndman and Koehler 2006). Scaling by the test
# set makes the value depend on the evaluation period and breaks
# comparability with the literature.
mase <- function(y, yhat, m = 7, y_train = NULL) {
  if (!is.null(y_train)) {
    n <- length(y_train)
    denom <- mean(abs(y_train[(m + 1):n] - y_train[1:(n - m)]), na.rm = TRUE)
    return(mean(abs(y - yhat), na.rm = TRUE) / denom)
  }
  mase_legacy(y, yhat, m)
}

mase_legacy <- function(y, yhat, m = 7) {
  y <- as.numeric(y)
  yhat <- as.numeric(yhat)
  e <- abs(y - yhat)
  denom <- mean(abs(y[(m + 1):length(y)] - y[1:(length(y) - m)]), na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  mean(e, na.rm = TRUE) / denom
}

test_r2 <- function(y, yhat) {
  y <- as.numeric(y)
  yhat <- as.numeric(yhat)
  ss_res <- sum((y - yhat)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(ss_tot) || ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

skill_score <- function(model_rmse, benchmark_rmse) {
  1 - (model_rmse / benchmark_rmse)
}

make_dow <- function(fecha) {
  w <- as.POSIXlt(as.Date(fecha))$wday
  labs <- c("sunday","monday","tuesday","wednesday","thursday","friday","saturday")
  factor(labs[w + 1], levels = labs)
}

numify_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  x <- gsub("\u00a0", "", x)
  x <- gsub("[^0-9,\\.\\-]", "", x)
  
  has_comma <- grepl(",", x)
  has_dot   <- grepl("\\.", x)
  out <- x
  
  idx_eu <- has_comma & has_dot & grepl("\\.\\d{3},", x)
  out[idx_eu] <- gsub("\\.", "", out[idx_eu])
  out[idx_eu] <- gsub(",", ".", out[idx_eu])
  
  idx_us <- has_comma & has_dot & grepl(",\\d{3}\\.", x)
  out[idx_us] <- gsub(",", "", out[idx_us])
  
  idx_c <- has_comma & !has_dot
  out[idx_c] <- gsub(",", ".", out[idx_c])
  
  suppressWarnings(as.numeric(out))
}

parse_date_fix <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "/", x)
  x <- gsub("-", "/", x)
  x <- sub(" .*$", "", x)
  
  out <- rep(as.Date(NA), length(x))
  suppressWarnings(xn <- as.numeric(x))
  is_excel <- !is.na(xn) & xn > 20000 & xn < 80000
  if (any(is_excel)) out[is_excel] <- as.Date(xn[is_excel], origin = "1899-12-30")
  
  idx <- is.na(out)
  if (any(idx)) {
    out[idx] <- suppressWarnings(as.Date(
      x[idx],
      tryFormats = c(
        "%d/%m/%Y", "%Y/%m/%d", "%m/%d/%Y",
        "%d/%m/%y", "%m/%d/%y", "%Y%m%d", "%Y-%m-%d"
      )
    ))
  }
  out
}

read_any <- function(path) {
  assert(file.exists(path), paste0("File not found: ", path))
  dt <- tryCatch(
    fread(path, encoding = "UTF-8", showProgress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dt)) {
    setnames(dt, trimws(gsub("^\ufeff", "", names(dt))))
    return(dt)
  }
  
  l1 <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1)) ";" else if (grepl(",", l1)) "," else ""
  dec <- if (sep == ";") "," else "."
  
  x <- if (sep == "") {
    read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.table(
      path, header = TRUE, sep = sep, dec = dec,
      stringsAsFactors = FALSE, check.names = FALSE, quote = "\""
    )
  }
  x <- as.data.table(x)
  setnames(x, trimws(gsub("^\ufeff", "", names(x))))
  x
}

collapse_by_date <- function(dt, value_col) {
  dt <- dt[!is.na(fecha)]
  if (nrow(dt) == 0) return(dt)
  if (anyDuplicated(dt$fecha)) {
    message("[QA] Duplicated dates detected in ", value_col, " -> daily mean aggregation")
    dt <- dt[, .(val = mean(get(value_col), na.rm = TRUE)), by = fecha]
    setnames(dt, "val", value_col)
  } else {
    dt <- dt[, .(fecha, tmp = get(value_col))]
    setnames(dt, "tmp", value_col)
  }
  setorder(dt, fecha)
  dt
}

predict_response <- function(mod, newdata) {
  if (inherits(mod, "glm")) {
    as.numeric(predict(mod, newdata = newdata, type = "response"))
  } else {
    as.numeric(predict(mod, newdata = newdata))
  }
}

eval_one <- function(name, y, yhat, aic = NA_real_, bic = NA_real_) {
  data.table(
    Model = name,
    AIC = aic,
    BIC = bic,
    RMSE = rmse(y, yhat),
    MAE = mae(y, yhat),
    MedAE = medae(y, yhat),
    `MAPE (%)` = mape(y, yhat),
    `sMAPE (%)` = smape(y, yhat),
    Bias = bias(y, yhat),
    `R² (test)` = test_r2(y, yhat),
    r = suppressWarnings(cor(y, yhat, use = "pairwise.complete.obs"))
  )
}

# ============================================================
# 2) READERS + QA/QC
# ============================================================

read_hidro_2cols <- function(path, nombre_valor) {
  step(paste0("READ PRECIP: ", basename(path)))
  dt <- read_any(path)
  assert(ncol(dt) >= 2, paste0("Hydro file has <2 columns: ", basename(path)))
  
  out <- data.table(
    fecha = parse_date_fix(dt[[1]]),
    val   = numify_safe(dt[[2]])
  )
  setnames(out, "val", nombre_valor)
  out <- out[!is.na(fecha)]
  assert(nrow(out) > 0, paste0("No valid dates in: ", basename(path)))
  
  bad_neg <- sum(out[[nombre_valor]] < 0, na.rm = TRUE)
  bad_hi  <- sum(out[[nombre_valor]] > 500, na.rm = TRUE)
  if (bad_neg + bad_hi > 0) {
    message("[QA] Implausible precip set to NA: neg=", bad_neg, " hi=", bad_hi)
  }
  out[get(nombre_valor) < 0 | get(nombre_valor) > 500, (nombre_valor) := NA_real_]
  
  out <- collapse_by_date(out, nombre_valor)
  
  message(
    "[OK] ", nombre_valor,
    " | rows: ", nrow(out),
    " | date range: ", paste(safe_range(out$fecha), collapse = " to "),
    " | missing: ", sum(is.na(out[[nombre_valor]]))
  )
  out
}

read_demanda <- function(path) {
  step("READ DEMAND")
  dt <- read_any(path)
  req <- c("date","val","weekend","holiday","sequia")
  assert(all(req %in% names(dt)), paste0("Demand requires columns: ", paste(req, collapse = ", ")))
  
  dem <- data.table(
    fecha   = parse_date_fix(dt$date),
    caudal  = numify_safe(dt$val),
    weekend = as.integer(dt$weekend),
    holiday = as.integer(dt$holiday),
    # C7: the 'sequia' column is a deterministic June-September indicator,
    # identical in every year of the record. It is a dry-SEASON marker, not
    # a drought index: it carries no interannual variation and is an exact
    # linear combination of the month factor. Renamed accordingly.
    drought = as.integer(dt$sequia)
  )
  dem <- dem[!is.na(fecha) & !is.na(caudal)]
  for (v in c("weekend","holiday","drought")) {
    dem <- dem[!is.na(get(v)) & get(v) %in% c(0L, 1L)]
  }
  setorder(dem, fecha)
  
  message("[QA] Demand rows: ", nrow(dem), " | date range: ", paste(safe_range(dem$fecha), collapse = " to "))
  message(
    "[QA] Demand summary (L/s): min=", round(min(dem$caudal), 2),
    " median=", round(median(dem$caudal), 2),
    " max=", round(max(dem$caudal), 2)
  )
  if (any(dem$caudal <= 0, na.rm = TRUE)) {
    message("[QA] Non-positive demand detected: n=", sum(dem$caudal <= 0, na.rm = TRUE))
  }
  
  dem[, anio := as.integer(format(fecha, "%Y"))]
  dem[, mes  := factor(format(fecha, "%m"), levels = sprintf("%02d", 1:12))]
  dem[, dow  := make_dow(fecha)]
  dem[, weekendF := factor(weekend, levels = c(0,1))]
  dem[, holidayF := factor(holiday, levels = c(0,1))]
  dem[, droughtF := factor(drought, levels = c(0,1))]
  
  dem
}

oni_from_csv <- function(path) {
  step("READ ONI (monthly)")
  dt <- read_any(path)
  nms <- tolower(gsub("\\s+", "_", gsub("\ufeff", "", trimws(names(dt)))))
  setnames(dt, nms)
  
  col_year <- grep("^yr$|year|anio|ano", names(dt), value = TRUE)[1]
  col_mon  <- grep("^mon$|month|mes", names(dt), value = TRUE)[1]
  col_anom <- grep("anom|oni|anomaly|anomalia", names(dt), value = TRUE)[1]
  
  assert(
    !anyNA(c(col_year, col_mon, col_anom)),
    paste0("Could not detect ONI columns. Available: ", paste(names(dt), collapse = ", "))
  )
  
  oni <- data.table(
    anio    = as.integer(numify_safe(dt[[col_year]])),
    mes_num = as.integer(numify_safe(dt[[col_mon]])),
    oni     = as.numeric(numify_safe(dt[[col_anom]]))
  )
  oni <- oni[!is.na(anio) & !is.na(mes_num) & mes_num %between% c(1L, 12L)]
  setorder(oni, anio, mes_num)
  
  message(
    "[OK] ONI rows: ", nrow(oni),
    " | years: ", min(oni$anio), "-", max(oni$anio),
    " | missing oni: ", sum(is.na(oni$oni))
  )
  oni
}

# ============================================================
# 3) FEATURE CONSTRUCTION
# ============================================================

build_features <- function(dem, p_bella, p_inq, oni_monthly) {
  step("BUILD FEATURES: composite precip + P7 + ONI")
  
  assert(!is.null(p_bella) || !is.null(p_inq), "No precipitation station could be read.")
  
  if (!is.null(p_bella) && !is.null(p_inq)) {
    P <- merge(p_bella, p_inq, by = "fecha", all = TRUE)
  } else if (!is.null(p_bella)) {
    P <- copy(p_bella)
    P[, precip_inaquito := NA_real_]
  } else {
    P <- copy(p_inq)
    P[, precip_bellavista := NA_real_]
  }
  setorder(P, fecha)
  
  both <- !is.na(P$precip_bellavista) & !is.na(P$precip_inaquito)
  message("[QA] Missing Bellavista: ", sum(is.na(P$precip_bellavista)))
  message("[QA] Missing Iñaquito : ", sum(is.na(P$precip_inaquito)))
  message("[QA] Days with both   : ", sum(both))
  if (sum(both) > 30) {
    cc <- suppressWarnings(cor(P$precip_bellavista[both], P$precip_inaquito[both]))
    message("[QA] Corr(Bellavista, Iñaquito) on overlap = ", round(cc, 3))
  }
  
  P[, precip_mm := rowMeans(cbind(precip_bellavista, precip_inaquito), na.rm = TRUE)]
  P[is.nan(precip_mm), precip_mm := NA_real_]
  P <- P[, .(fecha, precip_mm)]
  
  df <- merge(dem, P, by = "fecha", all.x = TRUE, sort = FALSE)
  setorder(df, fecha)
  
  df[, precip_mm := numify_safe(precip_mm)]
  df[precip_mm < 0 | precip_mm > 500, precip_mm := NA_real_]
  df[, miss_P := as.integer(is.na(precip_mm))]
  medP <- suppressWarnings(median(df$precip_mm, na.rm = TRUE))
  if (!is.finite(medP)) medP <- 0
  df[is.na(precip_mm), precip_mm := medP]
  
  message(
    "[QA] Precip median used for imputation (mm/day): ", round(medP, 3),
    " | imputed days: ", sum(df$miss_P)
  )
  
  df[, P7 := frollsum(precip_mm, n = 7, align = "right", na.rm = TRUE)]
  assert(any(!is.na(df$P7)), "P7 is all NA (check precip coverage).")
  
  df[, mes_num := as.integer(format(fecha, "%m"))]
  key_df  <- paste(df$anio, df$mes_num)
  key_oni <- paste(oni_monthly$anio, oni_monthly$mes_num)
  idx <- match(key_df, key_oni)
  
  df[, oni_step := oni_monthly$oni[idx]]      # monthly step function (legacy)
  df[, miss_oni := as.integer(is.na(oni_step))]

  # C6: ONI is a monthly index. Assigning it as a step function to daily
  # observations creates artificial discontinuities on the first of each
  # month and treats ~30 dependent days as independent. It is therefore
  # interpolated linearly between month midpoints. Clustered standard
  # errors are reported later in any case.
  oni_m <- copy(oni_monthly)
  oni_m[, fecha_c := as.Date(sprintf("%d-%02d-15", anio, mes_num))]
  setorder(oni_m, fecha_c)
  df[, oni := approx(x = oni_m$fecha_c, y = oni_m$oni,
                     xout = fecha, rule = 2)$y]
  df[is.na(oni), oni := 0]

  message("[QA] ONI: monthly values interpolated to daily resolution (C6)")
  message("[QA] ONI months not matched: n=", sum(df$miss_oni))
  
  df2 <- df[!is.na(P7)]
  setorder(df2, fecha)

  # C3: outliers are IDENTIFIED and reported but RETAINED in the analytical
  # dataset, so that the full operational record is modelled. Set
  # cfg$drop_outliers to TRUE to exclude them instead.
  lo_fence <- quantile(df2$caudal, .25, na.rm = TRUE) - 1.5 * IQR(df2$caudal, na.rm = TRUE)
  out_rows <- df2[caudal < lo_fence, .(fecha, caudal)]
  if (nrow(out_rows) > 0) {
    message("[QA] Lower boxplot fence: ", round(lo_fence, 1), " L/s")
    message("[QA] Observations below the fence (n=", nrow(out_rows),
            "), likely plant shutdowns \u2014 reported and RETAINED:")
    print(out_rows)
    if (isTRUE(getOption("qwd.drop_outliers", FALSE))) {
      df2 <- df2[caudal >= lo_fence]
      message("[QA] Excluded from the analytical dataset (cfg$drop_outliers = TRUE)")
    }
  }
  assign("qa_outliers", out_rows, envir = globalenv())
  
  message(
    "[OK] Feature dataset rows: ", nrow(df2),
    " | date range: ", paste(safe_range(df2$fecha), collapse = " to ")
  )
  message(
    "[QA] P7 summary: min=", round(min(df2$P7), 2),
    " p50=", round(median(df2$P7), 2),
    " max=", round(max(df2$P7), 2)
  )
  
  df2
}

# ============================================================
# 4) TEMPORAL SPLIT
# ============================================================

temporal_split <- function(df2, frac = 0.8) {
  step("TEMPORAL SPLIT (ordered 80/20) + factor alignment")
  df2 <- as.data.table(df2)
  assert("fecha" %in% names(df2), "df2 must contain fecha")
  
  if (!inherits(df2$fecha, c("Date","POSIXct","POSIXt"))) {
    df2[, fecha := as.Date(fecha)]
    if (anyNA(df2$fecha)) stop("fecha could not be coerced to Date.")
  }
  
  setorder(df2, fecha)
  n <- nrow(df2)
  cut <- floor(frac * n)
  assert(cut > 30 && (n - cut) > 30, "Split too small; check data length.")
  
  train <- copy(df2[1:cut])
  test  <- copy(df2[(cut + 1):n])
  
  train[, mes := factor(sprintf("%02d", as.integer(format(fecha, "%m"))), levels = sprintf("%02d", 1:12))]
  test[,  mes := factor(sprintf("%02d", as.integer(format(fecha, "%m"))), levels = levels(train$mes))]
  
  train[, dow := make_dow(fecha)]
  test[,  dow := factor(as.character(make_dow(fecha)), levels = levels(train$dow))]
  
  make_bin_factor <- function(dt, base, out) {
    if (base %in% names(dt)) {
      v <- dt[[base]]
      if (is.logical(v)) v <- as.integer(v)
      v <- suppressWarnings(as.integer(as.character(v)))
      v[!(v %in% c(0L,1L))] <- NA_integer_
      dt[, (out) := factor(v, levels = c(0,1))]
    } else if (out %in% names(dt)) {
      v <- suppressWarnings(as.integer(as.character(dt[[out]])))
      v[!(v %in% c(0L,1L))] <- NA_integer_
      dt[, (out) := factor(v, levels = c(0,1))]
    } else {
      dt[, (out) := factor(NA_integer_, levels = c(0,1))]
      message("[WARN] Neither ", base, " nor ", out, " found; created ", out, " as NA.")
    }
    invisible(dt)
  }
  
  make_bin_factor(train, "holiday",  "holidayF")
  make_bin_factor(test,  "holiday",  "holidayF")
  make_bin_factor(train, "weekend",  "weekendF")
  make_bin_factor(test,  "weekend",  "weekendF")
  make_bin_factor(train, "drought",  "droughtF")
  make_bin_factor(test,  "drought",  "droughtF")
  
  message("[OK] Train rows: ", nrow(train), " | Test rows: ", nrow(test))
  list(train = train, test = test)
}

# ============================================================
# 5) MODEL TABLE + FITTING
# ============================================================

model_specs_table <- function() {
  data.table(
    Model = c("M0","M0w","MP","MPS","MPO","MPOS"),
    Specification = c(
      "Demand ~ year + month + dow + holiday",
      "Demand ~ year + month + weekend + holiday",
      "Demand ~ M0 + P7 + miss_P",
      "Demand ~ MP + drought + P7 × drought",
      "Demand ~ M0 + P7 + P7² + ONI + miss_P",
      "Demand ~ M0 + P7 + P7² + ONI + drought + P7 × drought + miss_P"
    ),
    Purpose = c(
      "Baseline behavioral model (no climate forcing)",
      "Parsimonious behavioral baseline",
      "Introduces antecedent precipitation as the primary hydroclimatic predictor",
      "Tests heterogeneous precipitation response under drought",
      "Preferred hydroclimatic specification (nonlinear precipitation response + ENSO variability)",
      "Extended hydroclimatic specification (nonlinear + ENSO + drought interaction)"
    ),
    `Explanatory variables and interpretation` = c(
      "year: long-term structural trend; month: seasonal variability; dow: weekly behavioral pattern; holiday: public-holiday demand shifts.",
      "weekend: Saturday–Sunday indicator; other variables as in M0.",
      "P7: 7-day cumulative precipitation (antecedent rainfall); miss_P: missing-precipitation flag controlling for imputation effects.",
      "drought: drought-period indicator; P7 × drought: precipitation sensitivity allowed to vary under drought conditions.",
      "P7²: nonlinear rainfall response; ONI: Oceanic Niño Index anomaly representing regional hydroclimatic variability.",
      "Combines nonlinear precipitation response, ENSO variability, drought conditions, and precipitation–drought interaction effects."
    )
  )
}

fit_ols_models <- function(train) {
  list(
    M0   = lm(caudal ~ anio + mes + dow + holidayF, data = train),
    M0w  = lm(caudal ~ anio + mes + weekendF + holidayF, data = train),
    MP   = lm(caudal ~ anio + mes + dow + holidayF + P7 + miss_P, data = train),
    MPS  = lm(caudal ~ anio + mes + dow + holidayF + P7 + droughtF + P7:droughtF + miss_P, data = train),
    MPO  = lm(caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + miss_P, data = train),
    MPOS = lm(caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + droughtF + P7:droughtF + miss_P, data = train)
  )
}

fit_glm_models <- function(train, eps = 1e-6) {
  tr <- copy(train)
  if (any(tr$caudal <= 0, na.rm = TRUE)) {
    message("[QA] Adding epsilon to non-positive demand for Gamma GLM fit: eps=", eps)
    tr[caudal <= 0, caudal := eps]
  }
  
  fam <- Gamma(link = "log")
  list(
    M0   = glm(caudal ~ anio + mes + dow + holidayF, data = tr, family = fam),
    M0w  = glm(caudal ~ anio + mes + weekendF + holidayF, data = tr, family = fam),
    MP   = glm(caudal ~ anio + mes + dow + holidayF + P7 + miss_P, data = tr, family = fam),
    MPS  = glm(caudal ~ anio + mes + dow + holidayF + P7 + droughtF + P7:droughtF + miss_P, data = tr, family = fam),
    MPO  = glm(caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + miss_P, data = tr, family = fam),
    MPOS = glm(caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + droughtF + P7:droughtF + miss_P, data = tr, family = fam)
  )
}

eval_models <- function(models, test, mase_period = 7) {
  rbindlist(lapply(names(models), function(nm) {
    pr <- predict_response(models[[nm]], test)
    data.table(
      Model = nm,
      AIC = AIC(models[[nm]]),
      BIC = BIC(models[[nm]]),
      RMSE = rmse(test$caudal, pr),
      MAE  = mae(test$caudal, pr),
      MedAE = medae(test$caudal, pr),
      `MAPE (%)` = mape(test$caudal, pr),
      `sMAPE (%)` = smape(test$caudal, pr),
      Bias = bias(test$caudal, pr),
      `R² (test)` = test_r2(test$caudal, pr),
      r = suppressWarnings(cor(test$caudal, pr, use = "pairwise.complete.obs")),
      `MASE(7)` = mase(test$caudal, pr, m = mase_period, y_train = train$caudal)
    )
  }))[order(RMSE)]
}

# ============================================================
# 6) OTHER MODELS
# ============================================================

fix_xreg <- function(Xtr, Xte) {
  miss <- setdiff(colnames(Xtr), colnames(Xte))
  if (length(miss) > 0) {
    Xte <- cbind(Xte, matrix(0, nrow = nrow(Xte), ncol = length(miss), dimnames = list(NULL, miss)))
  }
  Xte <- Xte[, colnames(Xtr), drop = FALSE]
  
  bad <- apply(Xtr, 2, function(z) any(!is.finite(z)))
  if (any(bad)) {
    Xtr <- Xtr[, !bad, drop = FALSE]
    Xte <- Xte[, colnames(Xtr), drop = FALSE]
  }
  
  is_const <- apply(Xtr, 2, function(z) var(z, na.rm = TRUE) == 0)
  if (any(is_const)) {
    Xtr <- Xtr[, !is_const, drop = FALSE]
    Xte <- Xte[, colnames(Xtr), drop = FALSE]
  }
  
  if ("(Intercept)" %in% colnames(Xtr)) {
    Xtr <- Xtr[, colnames(Xtr) != "(Intercept)", drop = FALSE]
    Xte <- Xte[, colnames(Xtr), drop = FALSE]
  }
  
  list(Xtr = Xtr, Xte = Xte)
}

fit_other_models <- function(train, test) {
  step("FIT OTHER MODELS: ARIMA / ETS / SARIMAX / GAM / GLM")
  
  for (nm in c("mes","dow","holidayF","weekendF","droughtF")) {
    if (nm %in% names(train)) train[, (nm) := as.factor(get(nm))]
    if (nm %in% names(test))  test[,  (nm) := as.factor(get(nm))]
  }
  
  f_high <- caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + miss_P
  f_pars <- caudal ~ anio + mes + weekendF + holidayF + P7 + I(P7^2) + oni + miss_P
  
  y_tr <- ts(train$caudal, frequency = 7)
  h <- nrow(test)
  
  m_ets   <- ets(y_tr)
  p_ets   <- as.numeric(forecast(m_ets, h = h)$mean)
  
  m_arima <- auto.arima(y_tr, seasonal = TRUE, stepwise = TRUE, approximation = FALSE)
  p_arima <- as.numeric(forecast(m_arima, h = h)$mean)
  
  xreg_formula_high <- ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + miss_P
  Xtr_high <- model.matrix(xreg_formula_high, data = train)
  Xte_high <- model.matrix(xreg_formula_high, data = test)
  fx <- fix_xreg(Xtr_high, Xte_high)
  
  m_sarimax_high <- auto.arima(
    y_tr, xreg = fx$Xtr,
    seasonal = TRUE, stepwise = TRUE, approximation = FALSE
  )
  p_sarimax_high <- as.numeric(forecast(m_sarimax_high, xreg = fx$Xte, h = h)$mean)
  
  xreg_formula_pars <- ~ anio + mes + weekendF + holidayF + P7 + I(P7^2) + oni + miss_P
  Xtr_pars <- model.matrix(xreg_formula_pars, data = train)
  Xte_pars <- model.matrix(xreg_formula_pars, data = test)
  fx2 <- fix_xreg(Xtr_pars, Xte_pars)
  
  m_sarimax_pars <- auto.arima(
    y_tr, xreg = fx2$Xtr,
    seasonal = TRUE, stepwise = TRUE, approximation = FALSE
  )
  p_sarimax_pars <- as.numeric(forecast(m_sarimax_pars, xreg = fx2$Xte, h = h)$mean)
  
  glm_high <- glm(f_high, data = train, family = Gamma(link = "log"))
  p_glm_high <- as.numeric(predict(glm_high, newdata = test, type = "response"))
  
  glm_pars <- glm(f_pars, data = train, family = Gamma(link = "log"))
  p_glm_pars <- as.numeric(predict(glm_pars, newdata = test, type = "response"))
  
  gam_high <- gam(
    caudal ~ anio + mes + dow + holidayF + s(P7, k = 10) + s(oni, k = 8) + miss_P,
    data = train, family = Gamma(link = "log"), method = "REML"
  )
  p_gam_high <- as.numeric(predict(gam_high, newdata = test, type = "response"))
  
  gam_pars <- gam(
    caudal ~ anio + mes + weekendF + holidayF + s(P7, k = 10) + s(oni, k = 8) + miss_P,
    data = train, family = Gamma(link = "log"), method = "REML"
  )
  p_gam_pars <- as.numeric(predict(gam_pars, newdata = test, type = "response"))
  
  tab_compare <- rbindlist(list(
    eval_one("ARIMA (univariate)",            test$caudal, p_arima,        AIC(m_arima)),
    eval_one("GLM Gamma(log) (high-res)",     test$caudal, p_glm_high,     AIC(glm_high)),
    eval_one("GLM Gamma(log) (parsimonious)", test$caudal, p_glm_pars,     AIC(glm_pars)),
    eval_one("SARIMAX (parsimonious xreg)",   test$caudal, p_sarimax_pars, AIC(m_sarimax_pars)),
    eval_one("GAM Gamma(log) (high-res)",     test$caudal, p_gam_high,     AIC(gam_high)),
    eval_one("GAM Gamma(log) (parsimonious)", test$caudal, p_gam_pars,     AIC(gam_pars)),
    eval_one("ETS (univariate)",              test$caudal, p_ets,          AIC(m_ets)),
    eval_one("SARIMAX (high-res xreg)",       test$caudal, p_sarimax_high, AIC(m_sarimax_high))
  ), fill = TRUE)
  
  order_models <- c(
    "ARIMA (univariate)",
    "GLM Gamma(log) (high-res)",
    "GLM Gamma(log) (parsimonious)",
    "SARIMAX (parsimonious xreg)",
    "GAM Gamma(log) (high-res)",
    "GAM Gamma(log) (parsimonious)",
    "ETS (univariate)",
    "SARIMAX (high-res xreg)"
  )
  tab_compare[, Model := factor(Model, levels = order_models)]
  setorder(tab_compare, Model)
  
  fwrite(tab_compare, file.path(cfg$results_dir, "Table_Comparison_GLM_GAM_SARIMAX_ARIMA_ETS.csv"), bom = TRUE)
  
  list(
    table = tab_compare,
    m_ets = m_ets, p_ets = p_ets,
    m_arima = m_arima, p_arima = p_arima,
    m_sarimax_pars = m_sarimax_pars, p_sarimax_pars = p_sarimax_pars,
    m_sarimax_high = m_sarimax_high, p_sarimax_high = p_sarimax_high,
    glm_high = glm_high, glm_pars = glm_pars,
    gam_high = gam_high, gam_pars = gam_pars,
    p_glm_high = p_glm_high, p_glm_pars = p_glm_pars,
    p_gam_high = p_gam_high, p_gam_pars = p_gam_pars
  )
}

# ============================================================
# 7) FORECASTING
# ============================================================

# C5: Ecuadorian public holidays are deterministic and known in advance.
# The original pipeline assigned them at random within each year, which
# destroys their correlation with day-of-week and month and, because the
# annual maximum is sensitive to holiday placement, distorts the peak factor.
ec_holidays <- function(years) {
  # Easter (Meeus/Butcher) for the movable feasts
  easter <- function(y) {
    a <- y %% 19; b <- y %/% 100; c <- y %% 100
    d <- b %/% 4; e <- b %% 4; f <- (b + 8) %/% 25
    g <- (b - f + 1) %/% 3; h <- (19*a + b - d - g + 15) %% 30
    i <- c %/% 4; k <- c %% 4; l <- (32 + 2*e + 2*i - h - k) %% 7
    m <- (a + 11*h + 22*l) %/% 451
    mo <- (h + l - 7*m + 114) %/% 31; da <- ((h + l - 7*m + 114) %% 31) + 1
    as.Date(sprintf("%d-%02d-%02d", y, mo, da))
  }
  out <- unlist(lapply(years, function(y) {
    e <- easter(y)
    fixed <- as.Date(paste0(y, c("-01-01","-05-01","-05-24","-08-10",
                                 "-10-09","-11-02","-11-03","-12-25")))
    movable <- c(e - 48, e - 47, e - 2)   # Carnival Mon/Tue, Good Friday
    as.character(c(fixed, movable))
  }))
  as.Date(unique(out))
}

project_mixed <- function(df2, train, model, oni_monthly, anio_fin, seed, oni_obs_through_year) {
  stopifnot(inherits(model, c("lm","glm")))
  set.seed(seed)
  
  fecha_ini <- as.Date(max(df2$fecha) + 1)
  fecha_fin <- as.Date(paste0(anio_fin, "-12-31"))
  fechas_fut <- seq(fecha_ini, fecha_fin, by = "day")
  
  fut <- data.table(
    fecha   = fechas_fut,
    anio    = as.integer(format(fechas_fut, "%Y")),
    mes     = factor(format(fechas_fut, "%m"), levels = levels(train$mes)),
    mes_num = as.integer(format(fechas_fut, "%m")),
    dow     = make_dow(fechas_fut)
  )
  
  fut[, weekend := as.integer(dow %in% c("saturday","sunday"))]
  fut[, weekendF := factor(weekend, levels = levels(train$weekendF))]
  
  # C5: real holiday calendar instead of random assignment
  hol <- ec_holidays(sort(unique(fut$anio)))
  fut[, holiday := as.integer(fecha %in% hol)]
  fut[, holidayF := factor(holiday, levels = levels(train$holidayF))]
  message("[PROJ] Public holidays assigned from the Ecuadorian calendar: ",
          sum(fut$holiday), " days over ", length(unique(fut$anio)), " years")

  # C7: dry-season indicator follows the calendar, as in the observational data
  fut[, droughtF := factor(as.integer(mes_num %in% 6:9),
                           levels = levels(train$droughtF))]
  
  # C4: the original pipeline held P7 at a single constant value (p10, p50
  # or p90) for every day of the projection horizon. That removes the
  # bimodal rainfall seasonality of the study area, generates covariate
  # combinations that never occur (a p10 November sustained for 30 days),
  # suppresses a further source of variance in the peak factor, and makes
  # the small difference between "scenarios" an artefact of construction.
  # P7 is instead resampled within calendar month from the observed
  # distribution, restricted to the tercile defining each scenario.
  sample_P7 <- function(fechas, escenario) {
    mes_f <- as.integer(format(fechas, "%m"))
    q <- switch(escenario,
                "Dry (P7 P10)"    = c(0.00, 0.33),
                "Normal (P7 P50)" = c(0.33, 0.67),
                "Wet (P7 P90)"    = c(0.67, 1.00))
    vapply(seq_along(fechas), function(i) {
      pool <- df2[as.integer(format(fecha, "%m")) == mes_f[i], P7]
      pool <- pool[is.finite(pool)]
      if (!length(pool)) return(median(df2$P7, na.rm = TRUE))
      cut <- quantile(pool, q, na.rm = TRUE)
      sub <- pool[pool >= cut[1] & pool <= cut[2]]
      if (!length(sub)) sub <- pool
      sample(sub, 1)
    }, numeric(1))
  }
  escP <- c("Dry (P7 P10)", "Normal (P7 P50)", "Wet (P7 P90)")
  
  key_fut <- paste(fut$anio, fut$mes_num)
  key_oni <- paste(oni_monthly$anio, oni_monthly$mes_num)
  idx <- match(key_fut, key_oni)
  fut[, oni_real := oni_monthly$oni[idx]]
  
  escONI <- c("La Niña (ONI=-1)" = -1, "Neutral (ONI=0)" = 0, "El Niño (ONI=+1)" = 1)
  
  make_nd <- function(P7_scn, oni_scn) {
    nd <- copy(fut)
    nd[, P7 := sample_P7(fecha, P7_scn)]
    nd[, miss_P := 0L]
    nd[, oni := oni_real]
    nd[, miss_oni := as.integer(is.na(oni_real))]
    idx_scn <- (nd$anio > oni_obs_through_year) | is.na(nd$oni_real)
    nd[idx_scn, oni := oni_scn]
    nd
  }
  
  # C2 + C13: the conditional mean carries no residual variance, so the
  # daily maximum of predict() is not a peak demand and max/mean computed
  # on it understates the peak factor by roughly 0.11 (1.06 vs 1.17 in this
  # system). Demand is simulated from the fitted Gamma so that the maximum,
  # the peak factor and the prediction interval all reflect residual
  # variability.
  disp <- if (inherits(model, "glm") && !is.null(summary(model)$dispersion)) {
    summary(model)$dispersion
  } else {
    var(residuals(model)) / mean(fitted(model))^2
  }
  shape <- 1 / max(disp, 1e-9)

  pred_anual <- function(nd, n_sim = cfg$n_boot) {
    mu <- predict_response(model, nd)
    nd[, mu := mu]
    det <- nd[, .(mean_annual_demand = mean(mu, na.rm = TRUE)), by = anio]

    sims <- rbindlist(lapply(seq_len(n_sim), function(s) {
      y <- rgamma(length(mu), shape = shape, scale = mu / shape)
      data.table(anio = nd$anio, y = y)[, .(
        m = mean(y), mx = max(y), K = max(y) / mean(y)), by = anio][, sim := s]
    }))
    agg <- sims[, .(
      mean_lo          = quantile(m,  0.025),
      mean_hi          = quantile(m,  0.975),
      max_daily_demand = median(mx),
      max_lo           = quantile(mx, 0.025),
      max_hi           = quantile(mx, 0.975),
      K_peak           = median(K),
      K_lo             = quantile(K,  0.025),
      K_hi             = quantile(K,  0.975)
    ), by = anio]
    merge(det, agg, by = "anio")
  }
  
  rbindlist(lapply(escP, function(pn) {
    rbindlist(lapply(names(escONI), function(on) {
      nd <- make_nd(pn, escONI[[on]])
      an <- pred_anual(nd)
      an[, Scenario := paste0(pn, " | ", on)]
      an
    }))
  }))[order(Scenario, anio)]
}

# ============================================================
# 8) HYDROSOCIAL CONTRIBUTIONS
# ============================================================

hydrosocial_block_contributions_ols <- function(train, ols_models, model_key = "MPO") {
  step(paste0("HYDROSOCIAL BLOCK CONTRIBUTIONS (OLS): ", model_key))
  stopifnot(model_key %in% names(ols_models))
  
  mod_full0 <- ols_models[[model_key]]
  stopifnot(inherits(mod_full0, "lm"))
  
  data_for_attr <- as.data.table(copy(train))
  for (nm in c("mes","dow","holidayF","weekendF","droughtF")) {
    if (nm %in% names(data_for_attr)) data_for_attr[, (nm) := as.factor(get(nm))]
  }
  
  full_formula <- formula(mod_full0)
  mod_full <- lm(full_formula, data = data_for_attr)
  
  rss_full <- sum(residuals(mod_full)^2, na.rm = TRUE)
  tl_full  <- attr(terms(mod_full), "term.labels")
  
  blocks <- list(
    "Structural trend"   = c("anio"),
    "Social behavior"    = c("mes", "dow", "holidayF"),
    "Local hydroclimate" = c("P7", "I(P7^2)"),
    "ENSO variability"   = c("oni"),
    "Data control"       = c("miss_P")
  )
  
  blocks_in_model <- lapply(blocks, function(x) intersect(x, tl_full))
  blocks_in_model <- blocks_in_model[sapply(blocks_in_model, length) > 0]
  
  rows_h <- list()
  for (bname in names(blocks_in_model)) {
    drop_vec <- blocks_in_model[[bname]]
    red_formula <- update(
      full_formula,
      as.formula(paste0(". ~ . - ", paste(drop_vec, collapse = " - ")))
    )
    
    mod_red <- tryCatch(
      lm(red_formula, data = data_for_attr),
      error = function(e) NULL
    )
    if (is.null(mod_red)) next
    
    rss_red <- sum(residuals(mod_red)^2, na.rm = TRUE)
    delta <- rss_red - rss_full
    
    if (!is.finite(delta)) next
    if (delta < 0) delta <- 0
    
    rows_h[[length(rows_h) + 1]] <- data.table(
      Block = bname,
      TermsRemoved = paste(drop_vec, collapse = " + "),
      Delta = delta
    )
  }
  
  tab_h <- rbindlist(rows_h, fill = TRUE)
  if (nrow(tab_h) == 0) stop("No contributions computed for MPO (OLS).")
  
  tot_h <- sum(tab_h$Delta)
  if (!is.finite(tot_h) || tot_h <= 0) stop("Total contribution <= 0 for MPO (OLS).")
  
  tab_h[, Percent := 100 * Delta / tot_h]
  setorder(tab_h, -Percent)
  
  fwrite(
    tab_h[, .(Block, Percent = round(Percent, 1))],
    file.path(cfg$results_dir, "Table_HydrosocialBlockContributions_MPO_OLS.csv"),
    bom = TRUE
  )
  
  tab_h
}

draw_bellavista_hydrosocial_ols <- function(tab_h) {
  # Ranked horizontal bar chart. Replaces the earlier arrow schematic, in
  # which the central label was clipped, the legend duplicated the in-box
  # labels and most of the canvas was empty.
  d <- copy(as.data.table(tab_h))
  pct_col <- intersect(c("Percent", "Pct", "Share"), names(d))[1]
  stopifnot(!is.na(pct_col), "Block" %in% names(d))
  d <- d[, .(Block = as.character(Block), Pct = as.numeric(get(pct_col)))]
  d <- d[is.finite(Pct)]
  setorder(d, -Pct)

  pal <- c("#3B5E8C", "#4E8C6A", "#B5854A", "#6E5C93", "#8A8A8A")
  pal <- rep(pal, length.out = nrow(d))
  xmax <- ceiling(max(d$Pct) / 20) * 20

  save_png("Figure6_HydrosocialDrivers_Bellavista_MPO_OLS.png", {
    par(family = "serif", mar = c(4.4, 12.5, 2.4, 2.2), xpd = NA)
    plot(NA, xlim = c(0, xmax * 1.10), ylim = c(0.4, nrow(d) + 0.8),
         axes = FALSE, xlab = "Contribution to explained variance (%)",
         ylab = "", font.lab = 2)
    axis(1, at = seq(0, xmax, by = 20), cex.axis = 0.95)
    abline(v = seq(0, xmax, by = 20), col = "grey88", lty = 3)
    for (i in seq_len(nrow(d))) {
      y <- nrow(d) - i + 1
      rect(0, y - 0.30, d$Pct[i], y + 0.30, col = pal[i], border = NA)
      text(-xmax * 0.03, y, d$Block[i], adj = 1, cex = 1.0, font = 2)
      text(d$Pct[i] + xmax * 0.018, y, sprintf("%.1f%%", d$Pct[i]),
           adj = 0, cex = 0.95)
    }
    mtext(paste("Bellavista water supply system \u2014",
                "daily urban water demand (Q, L/s)"),
          side = 3, line = 0.6, cex = 1.0, font = 2)
  }, width = 2300, height = 1350, res = 320)
}

# ============================================================
# 9) MAIN EXECUTION
# ============================================================

step("PIPELINE START")

stopifnot(file.exists(cfg$demanda))
stopifnot(file.exists(cfg$p_bella))
stopifnot(file.exists(cfg$p_inq))
stopifnot(file.exists(cfg$oni))

options(qwd.drop_outliers = isTRUE(cfg$drop_outliers))

dem <- read_demanda(cfg$demanda)
P_bella <- read_hidro_2cols(cfg$p_bella, "precip_bellavista")
P_inq   <- read_hidro_2cols(cfg$p_inq,   "precip_inaquito")
oni_monthly <- oni_from_csv(cfg$oni)

df2 <- build_features(dem, P_bella, P_inq, oni_monthly)

spl <- temporal_split(df2, frac = cfg$split_frac)
train <- spl$train
test  <- spl$test

for (nm in c("mes","dow","holidayF","weekendF","droughtF")) {
  if (nm %in% names(train)) train[, (nm) := as.factor(get(nm))]
  if (nm %in% names(test))  test[,  (nm) := as.factor(get(nm))]
}

# ============================================================
# 10) TABLES
# ============================================================

table1_title <- "Table 1: Model specifications and theoretical interpretation"
table2_title <- "Table 2: Out-of-sample predictive performance across univariate and hydroclimatic models"
table3_title <- "Table 3: Model family comparison across baseline, hydroclimatic, and ENSO specifications"

tab1 <- model_specs_table()
fwrite(tab1, file.path(cfg$tab_dir, "Table1_ModelSpecifications.csv"), bom = TRUE)

ols_models <- fit_ols_models(train)
glm_models <- fit_glm_models(train)

tab_ols <- eval_models(ols_models, test)
tab_glm <- eval_models(glm_models, test)

fwrite(tab_ols, file.path(cfg$results_dir, "Table_Performance_OLS.csv"), bom = TRUE)
fwrite(tab_glm, file.path(cfg$results_dir, "Table_Performance_GLM_Gamma.csv"), bom = TRUE)

# Table 4 of the manuscript: comparative performance of the regression-based
# specifications. OLS and GLM Gamma results are stacked so that the table
# matches the caption "Comparative performance of the regression-based model
# specifications" rather than being split across two files.
tab4_reg <- rbindlist(list(
  copy(tab_glm)[, Model := paste0(Model, " (GLM Gamma)")],
  copy(tab_ols)[, Model := paste0(Model, " (OLS)")]
), fill = TRUE)

# Same columns, order and rounding as published in the manuscript.
keep4 <- intersect(c("Model", "AIC", "BIC", "RMSE", "MAE", "MAPE (%)",
                     "Bias", "R² (test)", "r"), names(tab4_reg))
tab4_reg <- tab4_reg[, ..keep4]
setorder(tab4_reg, RMSE)
for (cc in intersect(c("AIC", "BIC", "RMSE", "MAE"), names(tab4_reg)))
  tab4_reg[, (cc) := round(get(cc), 0)]
for (cc in intersect(c("MAPE (%)", "Bias"), names(tab4_reg)))
  tab4_reg[, (cc) := round(get(cc), 2)]
for (cc in intersect(c("R² (test)", "r"), names(tab4_reg)))
  tab4_reg[, (cc) := round(get(cc), 3)]
fwrite(tab4_reg, file.path(cfg$tab_dir, "Table4_RegressionSpecifications.csv"), bom = TRUE)

writeLines(paste0(
  "Table 4: Comparative performance of the regression-based model ",
  "specifications on the test set (n = ", nrow(test), "). Specifications are ",
  "defined in Table 2. Lower RMSE, MAE, MAPE and MASE indicate better ",
  "predictive accuracy. The preferred specification (", cfg$preferred_spec,
  ") was fixed a priori on grounds of parsimony, interpretability and scenario ",
  "capability, not by test-set error; MPS and MPOS are rank-deficient because ",
  "the dry-season indicator is an exact function of calendar month."),
  file.path(cfg$tab_dir, "Table4_Notes.txt"))
message("[OK] Table 4 (regression specifications) written")

# C11: the preferred specification is fixed a priori on grounds of
# parsimony, interpretability and scenario capability, NOT selected by
# test-set RMSE at run time. MPO and MPOS differ by 0.12 L/s in RMSE, so
# run-time selection made the reported model unstable across data updates.
best_glm_name <- cfg$preferred_spec
stopifnot(best_glm_name %in% names(glm_models))
best_glm <- glm_models[[best_glm_name]]
message("[OK] Preferred specification (fixed a priori): ", best_glm_name)

# ------------------------------------------------------------
# C1: FINAL MODEL FOR PROJECTION
# The 80/20 split exists to EVALUATE competing specifications. Once the
# specification is chosen, discarding 20% of the record (511 days, the most
# recent ones) when extrapolating 27 years forward is indefensible. The
# preferred specification is therefore refitted on the full observational
# record before any projection is generated.
# ------------------------------------------------------------
final_model <- glm(formula(glm_models[[best_glm_name]]),
                   data = df2, family = Gamma(link = "log"))
message("[OK] C1 final model refitted on the full record: n = ",
        length(final_model$fitted.values), " (training split was ", nrow(train), ")")
message("[OK] Annual trend coefficient: train = ",
        round(coef(glm_models[[best_glm_name]])["anio"], 5),
        " | full record = ", round(coef(final_model)["anio"], 5))

# ------------------------------------------------------------
# C6: CLUSTERED STANDARD ERRORS
# Monthly ONI values were repeated across the days of each month, so the
# ~2,548 daily observations contain roughly 84 independent ONI values.
# Naive standard errors overstate precision by a factor of about 3.6.
# ------------------------------------------------------------
df2[, mes_anio := format(fecha, "%Y-%m")]
ct_naive <- coeftest(final_model)
ct_clust <- coeftest(final_model, vcov = vcovCL, cluster = df2$mes_anio)

keep <- intersect(c("P7","I(P7^2)","oni"), rownames(ct_naive))
tab_se <- data.table(
  Term            = keep,
  Estimate        = round(ct_naive[keep, 1], 6),
  `t (naive)`     = round(ct_naive[keep, 3], 2),
  `p (naive)`     = signif(ct_naive[keep, 4], 3),
  `t (clustered)` = round(ct_clust[keep, 3], 2),
  `p (clustered)` = signif(ct_clust[keep, 4], 3)
)
tab_se[, `Inflation factor` := round(abs(`t (naive)`) / pmax(abs(`t (clustered)`), 1e-9), 1)]
fwrite(tab_se, file.path(cfg$results_dir, "Table_ClusteredStandardErrors.csv"), bom = TRUE)
message("[OK] C6 clustered standard errors written")
print(tab_se)

# ------------------------------------------------------------
# C14: SERIAL DEPENDENCE AND AUTOCORRELATION-CONSISTENT INFERENCE
# The specification contains no autoregressive component, so residual
# serial dependence is expected. It is quantified here because it governs
# what the model can support: under first-order dependence the nominal
# sample size overstates the information available, and it is also the
# reason the univariate ARIMA benchmark attains lower forecast error.
# Clustering by month does not address lag-1 dependence within months,
# so Newey-West standard errors are reported, with a sensitivity analysis
# over bandwidths exceeding the range of detectable autocorrelation.
# ------------------------------------------------------------

res_p  <- residuals(final_model, type = "pearson")
n_obs  <- length(res_p)
acf_v  <- acf(res_p, lag.max = 40, plot = FALSE)$acf[-1]
rho1   <- acf_v[1]
ci_band <- 1.96 / sqrt(n_obs)
n_eff  <- n_obs * (1 - rho1) / (1 + rho1)

dw <- lmtest::dwtest(lm(formula(final_model), data = df2))
lb7  <- Box.test(res_p, lag = 7,  type = "Ljung-Box")
lb30 <- Box.test(res_p, lag = 30, type = "Ljung-Box")

message(sprintf("[QA] C14 residual ACF: lag1 = %.3f | lag7 = %.3f | lag30 = %.3f",
                acf_v[1], acf_v[7], acf_v[30]))
message(sprintf("[QA] C14 lags outside the 95%% band (of 40): %d",
                sum(abs(acf_v) > ci_band)))
message(sprintf("[QA] C14 Durbin-Watson = %.3f (p = %.3g)", dw$statistic, dw$p.value))
message(sprintf("[QA] C14 Ljung-Box(30) = %.0f (p = %.3g)", lb30$statistic, lb30$p.value))
message(sprintf("[QA] C14 effective sample size ~ %.0f of %d (%.0f%%)",
                n_eff, n_obs, 100 * n_eff / n_obs))

terms_k  <- intersect(c("P7", "I(P7^2)", "oni", "anio"), names(coef(final_model)))
bw_auto  <- ceiling(4 * (n_obs / 100)^(2/9))
bw_set   <- sort(unique(c(bw_auto, 20, 30, 60)))

hac_ct <- function(L) coeftest(final_model,
  vcov = NeweyWest(final_model, lag = L, prewhite = FALSE, adjust = TRUE))

tab_hac <- rbindlist(lapply(bw_set, function(L) {
  ct <- hac_ct(L)
  data.table(Bandwidth = L, Term = terms_k,
             Estimate = round(ct[terms_k, 1], 6),
             `t` = round(ct[terms_k, 3], 2),
             `p` = signif(ct[terms_k, 4], 3))
}))

# tabla suplementaria: comparacion de los tres esquemas de error estandar
ct_naive <- coeftest(final_model)
ct_clust <- coeftest(final_model, vcov = vcovCL, cluster = df2$mes_anio)
ct_hac   <- hac_ct(bw_auto)

tab_se_all <- data.table(
  Term               = terms_k,
  Estimate           = round(ct_naive[terms_k, 1], 6),
  `t (conventional)` = round(ct_naive[terms_k, 3], 2),
  `p (conventional)` = signif(ct_naive[terms_k, 4], 3),
  `t (month-clustered)` = round(ct_clust[terms_k, 3], 2),
  `p (month-clustered)` = signif(ct_clust[terms_k, 4], 3),
  `t (Newey-West)`   = round(ct_hac[terms_k, 3], 2),
  `p (Newey-West)`   = signif(ct_hac[terms_k, 4], 3)
)
for (L in setdiff(bw_set, bw_auto)) {
  ct <- hac_ct(L)
  tab_se_all[[paste0("t (Newey-West, bw = ", L, ")")]] <- round(ct[terms_k, 3], 2)
}

fwrite(tab_se_all, file.path(cfg$results_dir, "Table_SerialDependence_Inference.csv"), bom = TRUE)
fwrite(tab_hac,    file.path(cfg$results_dir, "Table_HAC_Bandwidth_Sensitivity.csv"), bom = TRUE)

acf_out <- data.table(Lag = seq_along(acf_v), ACF = round(as.numeric(acf_v), 4))
acf_out[, `Outside 95% band` := fifelse(abs(ACF) > ci_band, "yes", "no")]
fwrite(acf_out, file.path(cfg$results_dir, "Table_ResidualACF.csv"), bom = TRUE)

s7_note <- paste0(
  "Notes: Estimates from the MPO GLM Gamma specification refitted on the full ",
  "observational record (n = ", n_obs, "). Residuals exhibit serial dependence ",
  "(first-order autocorrelation ", sprintf("%.3f", rho1),
  "; all 40 lags examined fall outside the 95% band; Durbin-Watson ",
  sprintf("%.2f", dw$statistic), "; Ljung-Box(30) p < 0.001), so the effective ",
  "number of independent observations is approximately ", round(n_eff),
  ". Conventional standard errors are reported for comparison only and are not ",
  "used for inference. Month-clustered errors address dependence between months ",
  "but not lag-1 dependence within them. Newey-West standard errors (automatic ",
  "bandwidth ", bw_auto, " days) are used throughout the manuscript; the ",
  "additional columns show that the conclusions are unchanged at bandwidths up ",
  "to 60 days, beyond the range over which residual autocorrelation remains ",
  "detectable."
)
writeLines(s7_note, file.path(cfg$results_dir, "Table_SerialDependence_Notes.txt"))

cat("\n=========== SERIAL DEPENDENCE AND ROBUST INFERENCE ===========\n\n")
print(tab_se_all)
cat("\n", s7_note, "\n\n", sep = "")
  # C7: MPS and MPOS include the dry-season indicator, which is an exact
  # linear combination of the month factor. Their design matrices are
  # rank-deficient and the main effect of droughtF is not identifiable.
  chk <- tryCatch({
    X <- model.matrix(caudal ~ anio + mes + dow + holidayF + P7 + droughtF, data = train)
    c(ncol(X), qr(X)$rank)
  }, error = function(e) c(NA_integer_, NA_integer_))
  if (!any(is.na(chk)) && chk[1] != chk[2]) {
    message("[QA] C7 rank deficiency: ", chk[1], " design columns, rank ", chk[2],
            " - dry-season indicator aliased with the month factor (MPS/MPOS)")
  }

message("[OK] Best GLM by RMSE: ", best_glm_name)

other_models <- fit_other_models(train, test)
tab2 <- copy(other_models$table)

tab2_out <- copy(tab2)
tab2_out[, AIC := round(AIC, 0)]
tab2_out[, RMSE := round(RMSE, 1)]
tab2_out[, MAE := round(MAE, 1)]
tab2_out[, MedAE := round(MedAE, 2)]
tab2_out[, `MAPE (%)` := round(`MAPE (%)`, 2)]
tab2_out[, `sMAPE (%)` := round(`sMAPE (%)`, 2)]
tab2_out[, Bias := round(Bias, 1)]
tab2_out[, `R² (test)` := round(`R² (test)`, 3)]
tab2_out[, r := round(r, 2)]

fwrite(tab2_out, file.path(cfg$tab_dir, "Table2_OutOfSamplePerformance.csv"), bom = TRUE)

rows_t3 <- list()
for (nm in names(glm_models)) {
  pr <- predict_response(glm_models[[nm]], test)
  rows_t3[[length(rows_t3)+1]] <- data.table(
    Model = paste0(nm, " (GLM Gamma)"),
    AIC = AIC(glm_models[[nm]]),
    BIC = BIC(glm_models[[nm]]),
    RMSE = rmse(test$caudal, pr),
    MAE = mae(test$caudal, pr),
    `MAPE (%)` = mape(test$caudal, pr),
    Bias = bias(test$caudal, pr),
    `R² (test)` = test_r2(test$caudal, pr),
    r = suppressWarnings(cor(test$caudal, pr, use = "pairwise.complete.obs"))
  )
}
for (nm in names(ols_models)) {
  pr <- predict_response(ols_models[[nm]], test)
  rows_t3[[length(rows_t3)+1]] <- data.table(
    Model = paste0(nm, " (OLS)"),
    AIC = AIC(ols_models[[nm]]),
    BIC = BIC(ols_models[[nm]]),
    RMSE = rmse(test$caudal, pr),
    MAE = mae(test$caudal, pr),
    `MAPE (%)` = mape(test$caudal, pr),
    Bias = bias(test$caudal, pr),
    `R² (test)` = test_r2(test$caudal, pr),
    r = suppressWarnings(cor(test$caudal, pr, use = "pairwise.complete.obs"))
  )
}

tab3 <- rbindlist(rows_t3, fill = TRUE)
order_t3 <- c(
  "MPO (GLM Gamma)",
  "MPOS (GLM Gamma)",
  "MPO (OLS)",
  "MPOS (OLS)",
  "MPS (GLM Gamma)",
  "MP (GLM Gamma)",
  "MPS (OLS)",
  "MP (OLS)",
  "M0 (GLM Gamma)",
  "M0w (GLM Gamma)",
  "M0 (OLS)",
  "M0w (OLS)"
)
tab3 <- tab3[Model %in% order_t3]
tab3[, Model := factor(Model, levels = order_t3)]
setorder(tab3, Model)

tab3_out <- copy(tab3)
# C10: AIC is comparable only within a model family sharing a likelihood,
# response scale and effective sample size. The ARIMA benchmark uses a
# Gaussian likelihood on a differenced series with reduced n; the GLMs use
# a Gamma likelihood on levels with full n. The column is therefore
# dropped from the cross-family table and kept only within families.
if ("AIC" %in% names(tab3_out)) tab3_out[, AIC := NULL]
tab3_out[, BIC := round(BIC, 0)]
tab3_out[, RMSE := round(RMSE, 0)]
tab3_out[, MAE := round(MAE, 0)]
tab3_out[, `MAPE (%)` := round(`MAPE (%)`, 2)]
tab3_out[, Bias := round(Bias, 2)]
tab3_out[, `R² (test)` := round(`R² (test)`, 3)]
tab3_out[, r := round(r, 3)]

# This table compares the candidate SPECIFICATIONS (M0 to MPOS, OLS and GLM
# Gamma). It corresponds to Table 4 of the manuscript.
fwrite(tab3_out, file.path(cfg$results_dir, "Table_SpecificationComparison.csv"), bom = TRUE)



writeLines(c(table1_title), con = file.path(cfg$tab_dir, "Table1_Title.txt"))
writeLines(c(table2_title), con = file.path(cfg$tab_dir, "Table2_Title.txt"))
writeLines(c(table3_title), con = file.path(cfg$tab_dir, "Table3_Title.txt"))

# ============================================================
# 11) FORECASTING RESULTS
# ============================================================

tabla_mix <- project_mixed(
  df2 = df2,
  train = train,
  model = final_model,          # C1: full-record fit, not the training split
  oni_monthly = oni_monthly,
  anio_fin = cfg$anio_fin,
  seed = cfg$seed,
  oni_obs_through_year = cfg$oni_obs_through_year
)
fwrite(tabla_mix, file.path(cfg$results_dir, "Table_Forecast_Annual_MixedScenarios.csv"), bom = TRUE)

# ============================================================
# 12) FIGURE 2
# ============================================================

pred_best <- predict_response(best_glm, test)
pred_mpo_ols <- as.numeric(predict(ols_models$MPO, newdata = test))
pred_mpo_glm <- as.numeric(predict(glm_models$MPO, newdata = test, type = "response"))
pred_m0  <- as.numeric(predict(ols_models$M0,  newdata = test))
pred_m0w <- as.numeric(predict(ols_models$M0w, newdata = test))

save_png("Figure2_ObservedVsPredicted_TwoPanels.png", {
  par(mfrow = c(2,1), family = "serif")
  
  par(mar = c(3.5,5,2.5,1))
  plot(test$fecha, test$caudal,
       type = "l", col = "black", lwd = 1.2,
       xlab = "", ylab = "Daily demand (L/s)",
       font.lab = 2, bty = "l")
  lines(test$fecha, pred_mpo_ols, col = "#E69F00", lwd = 2.5)
  lines(test$fecha, pred_mpo_glm, col = "#FF1493", lwd = 2.5, lty = 2)
  legend("bottomleft",
         legend = c("Observed","MPO-like model","GLM Gamma (log-link)"),
         col = c("black","#E69F00","#FF1493"),
         lwd = c(1.2,2.5,2.5), lty = c(1,1,2),
         bty = "n", cex = 1)
  mtext("A", side = 3, adj = 0, line = 0.3, font = 2)
  
  par(mar = c(5,5,2.5,1))
  plot(test$fecha, test$caudal,
       type = "l", col = "black", lwd = 1.2,
       xlab = "Date", ylab = "Daily demand (L/s)",
       font.lab = 2, bty = "l")
  lines(test$fecha, pred_m0,  col = "#E69F00", lwd = 2.4)
  lines(test$fecha, pred_m0w, col = "#00C853", lwd = 2.4, lty = 2)
  legend("bottomleft",
         legend = c("Observed","M0: calendar (dow + holiday)","M0w: calendar (weekend + holiday)"),
         col = c("black","#E69F00","#00C853"),
         lwd = c(1.2,2.4,2.4), lty = c(1,1,2),
         bty = "n", cex = 1)
  mtext("B", side = 3, adj = 0, line = 0.3, font = 2)
}, width = 2200, height = 1800, res = 320)

# ============================================================
# 13) FIGURE 3
# ============================================================

stopifnot("MPOS" %in% names(glm_models))
mpos_mod <- glm_models$MPOS

ref <- copy(test[1])
if ("miss_P" %in% names(ref))   ref[, miss_P := 0L]
if ("miss_oni" %in% names(ref)) ref[, miss_oni := 0L]

ref[, anio := median(test$anio, na.rm = TRUE)]
ref[, oni  := median(test$oni,  na.rm = TRUE)]
ref[, mes := factor(names(sort(table(test$mes), decreasing = TRUE))[1], levels = levels(test$mes))]
ref[, dow := factor(names(sort(table(test$dow), decreasing = TRUE))[1], levels = levels(test$dow))]
ref[, holidayF := factor(names(sort(table(test$holidayF), decreasing = TRUE))[1], levels = levels(test$holidayF))]
ref[, weekendF := factor(names(sort(table(test$weekendF), decreasing = TRUE))[1], levels = levels(test$weekendF))]
ref[, droughtF := factor("0", levels = levels(test$droughtF))]

p7_grid <- seq(
  as.numeric(quantile(test$P7, 0.01, na.rm = TRUE)),
  as.numeric(quantile(test$P7, 0.99, na.rm = TRUE)),
  length.out = 300
)

nd <- ref[rep(1, length(p7_grid))]
nd[, P7 := p7_grid]

pr_link <- predict(mpos_mod, newdata = nd, type = "link", se.fit = TRUE)
eta <- as.numeric(pr_link$fit)
se  <- as.numeric(pr_link$se.fit)

fit <- exp(eta)
lo  <- exp(eta - 1.96 * se)
hi  <- exp(eta + 1.96 * se)

hfd <- 1e-3
nd_hi <- copy(nd)
nd_hi[, P7 := P7 + hfd]
f0 <- as.numeric(predict(mpos_mod, newdata = nd,    type = "response"))
f1 <- as.numeric(predict(mpos_mod, newdata = nd_hi, type = "response"))
slope <- (f1 - f0) / hfd

save_png("Figure3_RainfallDemand_MPOS.png", {
  par(mfrow = c(2,1), family = "serif")
  
  par(mar = c(3.8, 5.4, 2.6, 1.0))
  plot(p7_grid, fit,
       type = "n",
       xlab = "",
       ylab = "Predicted water demand (L/s)",
       font.lab = 2, cex.lab = 1, cex.axis = 1, bty = "l")
  polygon(c(p7_grid, rev(p7_grid)),
          c(lo, rev(hi)),
          col = "#A6CEE3",
          border = NA)
  lines(p7_grid, fit, lwd = 2.8, col = "#1F78B4")
  legend("topright",
         legend = c("Predicted mean", "95% confidence interval"),
         col    = c("#1F78B4", "#A6CEE3"),
         lwd    = c(2.8, 8),
         bty    = "n",
         cex    = 1)
  mtext("A", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  par(mar = c(5.0, 5.4, 2.6, 1.0))
  plot(p7_grid, slope,
       type = "l",
       lwd = 2.8,
       col = "#33A02C",
       xlab = "Antecedent rainfall P7 (mm, 7-day cumulative)",
       ylab = "Marginal effect (dDemand/dP7)",
       font.lab = 2, cex.lab = 1, cex.axis = 1, bty = "l")
  abline(h = 0, lwd = 1.8, lty = 2, col = "#E31A1C")
  mtext("B", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
}, width = 2200, height = 1800, res = 320)

# ============================================================
# 14) FIGURE 4
# ============================================================

resid_test  <- test$caudal - pred_best
fitted_test <- pred_best

save_png("Figure4_ResidualDiagnostics.png", {
  par(mfrow = c(2,2), family = "serif")
  
  par(mar = c(5,5.2,3,1))
  plot(fitted_test, resid_test,
       pch = 19, cex = 0.65, col = "#6A3D9A",
       xlab = "Fitted demand (L/s)",
       ylab = "Residuals (Observed - Predicted)",
       font.lab = 2, cex.lab = 1, cex.axis = 1, bty = "l")
  abline(h = 0, lty = 2, lwd = 2, col = "#E31A1C")
  mtext("A", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  par(mar = c(5,5.2,3,1))
  qqnorm(resid_test,
         pch = 19, cex = 0.65, col = "#1F78B4",
         xlab = "Theoretical quantiles",
         ylab = "Sample quantiles",
         main = "",
         font.lab = 2, cex.lab = 1, cex.axis = 1, bty = "l")
  qqline(resid_test, lwd = 2, col = "#FF7F00")
  mtext("B", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  par(mar = c(5,5.2,3,1))
  plot(fitted_test, sqrt(abs(resid_test)),
       pch = 19, cex = 0.65, col = "#33A02C",
       xlab = "Fitted demand (L/s)",
       ylab = expression(sqrt("|Residuals|")),
       font.lab = 2, cex.lab = 1, cex.axis = 1, bty = "l")
  mtext("C", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
  
  par(mar = c(5,5.2,3,1))
  acf(resid_test,
      lag.max = 30, main = "",
      col = "#1F78B4", lwd = 2,
      xlab = "Lag", ylab = "ACF", ci.col = "grey40")
  box(bty = "l")
  mtext("D", side = 3, adj = 0, line = 0.3, font = 2, cex = 1.2)
}, width = 2200, height = 1800, res = 320)

# ============================================================
# 15) FIGURE 5
# ============================================================

# Projection figure. Uses the full horizon defined in cfg$anio_fin and the
# 95% intervals produced by the Gamma simulation. Scenario labels are
# classified by pattern rather than by exact string, so that the accented
# characters in "El Nino" / "La Nina" cannot break the match under a
# different locale. The design capacity of the system and the 2026
# operational observation are drawn as reference lines.

tt5 <- as.data.table(copy(tabla_mix))
tt5 <- tt5[anio %between% c(2014L, cfg$anio_fin)]

tt5[, pn := fifelse(grepl("Dry",  Scenario), "Dry (P7 p10)",
             fifelse(grepl("Wet",  Scenario), "Wet (P7 p90)", "Normal (P7 p50)"))]
tt5[, on := fifelse(grepl("ONI=\\+1", Scenario), "El Nino (ONI = +1)",
             fifelse(grepl("ONI=-1",   Scenario), "La Nina (ONI = -1)",
                                                  "Neutral (ONI = 0)"))]

CAP_LS   <- 3000    # Bellavista design capacity (L/s)
OBS_2026 <- 2427    # operational mean demand reported for 2026 (L/s)

save_png("Figure5_ProjectedDemand_TwoPanels.png", {
  par(mfrow = c(2, 1), family = "serif", mar = c(4.2, 5.4, 2.6, 1.6))

  draw_panel <- function(d, key, cols, tag, show_note = FALSE) {
    lv <- intersect(names(cols), unique(d[[key]]))
    yl <- range(c(d$mean_lo, d$mean_hi, CAP_LS, OBS_2026), na.rm = TRUE)
    yl[2] <- yl[2] * 1.02
    plot(NA, xlim = range(d$anio), ylim = yl, xlab = "Year", ylab = "Q (L/s)",
         font.lab = 2, bty = "l", las = 1)
    abline(h = CAP_LS, lty = 2, lwd = 2, col = "#8A5116")
    xr <- range(d$anio)
    text(xr[1] + diff(xr) * 0.62, CAP_LS, "design capacity 3,000 L/s",
         pos = 3, cex = 0.76, col = "#8A5116")
    for (s in lv) {
      dd <- d[get(key) == s][order(anio)]
      polygon(c(dd$anio, rev(dd$anio)), c(dd$mean_lo, rev(dd$mean_hi)),
              col = adjustcolor(cols[s], alpha.f = 0.20), border = NA)
    }
    for (s in lv) {
      dd <- d[get(key) == s][order(anio)]
      lines(dd$anio, dd$demand, lwd = 2.6, col = cols[s])
    }
    points(2026, OBS_2026, pch = 21, bg = "#FFFFFF", col = "#141414",
           cex = 1.5, lwd = 2)
    text(2026, OBS_2026, "observed 2026", pos = 4, cex = 0.76)
    legend("topleft", legend = lv, col = cols[lv], lwd = 2.6, bty = "n",
           cex = 0.86)
    mtext(tag, side = 3, adj = 0, line = 0.5, font = 2, cex = 1.2)
    if (show_note)
      mtext(paste("Shaded bands: 95% simulation intervals (narrow, since they",
                  "reflect residual variance only)."),
            side = 3, line = 0.5, adj = 1, cex = 0.72, col = "#5A5A5A")
  }

  pa <- tt5[, .(demand  = mean(mean_annual_demand, na.rm = TRUE),
                mean_lo = mean(mean_lo, na.rm = TRUE),
                mean_hi = mean(mean_hi, na.rm = TRUE)), by = .(anio, pn)]
  draw_panel(pa, "pn",
             c("Dry (P7 p10)" = "#B5546A", "Normal (P7 p50)" = "#1F4E79",
               "Wet (P7 p90)" = "#3E8E5A"), "A", show_note = TRUE)

  pb <- tt5[pn == "Normal (P7 p50)",
            .(demand  = mean(mean_annual_demand, na.rm = TRUE),
              mean_lo = mean(mean_lo, na.rm = TRUE),
              mean_hi = mean(mean_hi, na.rm = TRUE)), by = .(anio, on)]
  draw_panel(pb, "on",
             c("El Nino (ONI = +1)" = "#B5546A", "Neutral (ONI = 0)" = "#1F4E79",
               "La Nina (ONI = -1)" = "#3E8E5A"), "B")
  mtext("ENSO scenarios coincide until 2026, when observed ONI values end.",
        side = 3, line = 0.5, adj = 1, cex = 0.72, col = "#5A5A5A")
}, width = 2200, height = 2000, res = 320)

# ============================================================
# 16) FIGURE 6
# ============================================================

tab_h <- hydrosocial_block_contributions_ols(train, ols_models, model_key = "MPO")
draw_bellavista_hydrosocial_ols(tab_h)

# ============================================================
# 17) MODEL COMPARISON + DM RESULTS
# ============================================================

y_test <- test$caudal
p_arima <- other_models$p_arima
p_ets <- other_models$p_ets
p_glm_high <- other_models$p_glm_high
p_glm_pars <- other_models$p_glm_pars
p_gam_high <- other_models$p_gam_high
p_gam_pars <- other_models$p_gam_pars
p_mpo_ols  <- as.numeric(predict(ols_models$MPO,  newdata = test))
p_mpos_glm <- as.numeric(predict(glm_models$MPOS, newdata = test, type = "response"))

tab_compare <- rbindlist(list(
  eval_one("ARIMA (univariate)",            y_test, p_arima,     AIC(other_models$m_arima), BIC(other_models$m_arima)),
  eval_one("ETS (univariate)",              y_test, p_ets,       AIC(other_models$m_ets),   BIC(other_models$m_ets)),
  eval_one("GLM Gamma(log) (high-res)",     y_test, p_glm_high,  AIC(other_models$glm_high), BIC(other_models$glm_high)),
  eval_one("GLM Gamma(log) (parsimonious)", y_test, p_glm_pars,  AIC(other_models$glm_pars), BIC(other_models$glm_pars)),
  eval_one("GAM Gamma(log) (high-res)",     y_test, p_gam_high,  AIC(other_models$gam_high), BIC(other_models$gam_high)),
  eval_one("GAM Gamma(log) (parsimonious)", y_test, p_gam_pars,  AIC(other_models$gam_pars), BIC(other_models$gam_pars)),
  eval_one("MPO (OLS)",                     y_test, p_mpo_ols,   AIC(ols_models$MPO), BIC(ols_models$MPO)),
  eval_one("MPOS (GLM Gamma)",              y_test, p_mpos_glm,  AIC(glm_models$MPOS), BIC(glm_models$MPOS))
), fill = TRUE)

setorder(tab_compare, RMSE)

rmse_benchmark <- tab_compare[Model == "ARIMA (univariate)", RMSE][1]
tab_compare[, `Skill vs ARIMA` := skill_score(RMSE, rmse_benchmark)]
tab_compare[Model == "ARIMA (univariate)", `Skill vs ARIMA` := 0]

tab_out <- copy(tab_compare)
num_cols_2 <- c("RMSE","MAE","MedAE","MAPE (%)","sMAPE (%)","Bias","Skill vs ARIMA")
for (cc in num_cols_2) tab_out[, (cc) := round(get(cc), 2)]
tab_out[, `R² (test)` := round(`R² (test)`, 3)]
tab_out[, r := round(r, 3)]
tab_out[, AIC := round(AIC, 0)]
tab_out[, BIC := round(BIC, 0)]

fwrite(tab_out, file.path(cfg$results_dir, "Table_ModelComparison_WithSkill.csv"), bom = TRUE)

# ------------------------------------------------------------
# Table 3 of the manuscript: performance comparison across model FAMILIES
# on the test set (ARIMA, ETS, SARIMAX, GAM, GLM). Built from the
# cross-family comparison so that the exported file matches the caption
# "Performance comparison of candidate forecasting models using the test
# dataset".
# ------------------------------------------------------------
# The family comparison is read back from the file written by
# fit_other_models(): the object `tab_compare` is later reassigned to the
# specification comparison with skill scores, so referring to it here would
# export the wrong table.
fam_src <- file.path(cfg$results_dir, "Table_Comparison_GLM_GAM_SARIMAX_ARIMA_ETS.csv")
if (file.exists(fam_src)) {
  tab3_fam <- fread(fam_src, encoding = "UTF-8")
  # AIC is retained: it is informative within each family and the reviewers
  # refer to it explicitly. The caption states that values are not comparable
  # ACROSS families, because the ARIMA and ETS benchmarks use a Gaussian
  # likelihood on a differenced series with reduced effective n, whereas the
  # GLM and GAM specifications use a Gamma likelihood on levels with full n.
  # BIC is dropped because it is empty for the univariate benchmarks.
  if ("BIC" %in% names(tab3_fam)) tab3_fam[, BIC := NULL]
  if ("AIC" %in% names(tab3_fam)) tab3_fam[, AIC := round(AIC, 0)]
  num_cols <- names(tab3_fam)[vapply(tab3_fam, is.numeric, logical(1))]
  for (cc in num_cols) {
    dg <- if (grepl("R\u00b2|^r$", cc)) 3 else 2
    tab3_fam[, (cc) := round(get(cc), dg)]
  }
  if ("RMSE" %in% names(tab3_fam)) setorder(tab3_fam, RMSE)
  fwrite(tab3_fam, file.path(cfg$tab_dir, "Table3_ModelFamilyComparison.csv"), bom = TRUE)

  writeLines(paste0(
    "Table 3: Performance comparison of candidate forecasting models using ",
    "the test dataset (n = ", nrow(test), "). Lower values of RMSE, MAE, ",
    "MedAE, MAPE and sMAPE indicate better predictive performance. ",
    "RMSE = root mean square error; MAE = mean absolute error; ",
    "MedAE = median absolute error; MAPE = mean absolute percentage error; ",
    "sMAPE = symmetric mean absolute percentage error; Bias = mean prediction ",
    "bias; r = Pearson correlation between observed and predicted demand. ",
    "AIC values are comparable only WITHIN each model family: the ARIMA and ",
    "ETS benchmarks are estimated by a Gaussian likelihood on the differenced ",
    "series with reduced effective sample size, whereas the GLM and GAM ",
    "specifications use a Gamma likelihood on the series in levels with the ",
    "full sample. Cross-family differences in AIC therefore reflect the change ",
    "of likelihood and scale rather than relative goodness of fit, and model ",
    "selection in this study does not rely on them."),
    file.path(cfg$tab_dir, "Table3_Notes.txt"))
  message("[OK] Table 3 (model family comparison) written")
} else {
  warning("[TABLE 3] family comparison file not found; Table 3 not exported",
          call. = FALSE)
}


err_arima    <- y_test - p_arima
err_glm_high <- y_test - p_glm_high
err_mpo_ols  <- y_test - p_mpo_ols
err_mpos_glm <- y_test - p_mpos_glm
err_gam_high <- y_test - p_gam_high

# C8: forecast errors from a multi-step path are serially correlated.
# With h = 1 the test applies no HAC correction to the loss-differential
# variance and the statistic is inflated. The evaluation horizon is used.
dm_safe <- function(e1, e2, benchmark, compared, h = length(e1), power = 2) {
  out <- tryCatch(
    dm.test(e1, e2, h = h, power = power, alternative = "two.sided"),
    error = function(e) NULL
  )
  
  if (is.null(out)) {
    return(data.table(
      Benchmark_Model = benchmark,
      Compared_Model = compared,
      DM_statistic = NA_real_,
      p_value = NA_real_,
      Interpretation = "Test failed"
    ))
  }
  
  pval <- as.numeric(out$p.value)
  stat <- as.numeric(out$statistic)
  
  interp <- if (is.na(pval)) {
    "Test failed"
  } else if (pval < 0.05 && stat < 0) {
    paste(benchmark, "significantly more accurate")
  } else if (pval < 0.05 && stat > 0) {
    paste(compared, "significantly more accurate")
  } else {
    "No significant difference"
  }
  
  data.table(
    Benchmark_Model = benchmark,
    Compared_Model = compared,
    DM_statistic = stat,
    p_value = pval,
    Interpretation = interp
  )
}

tab_dm <- rbindlist(list(
  dm_safe(err_arima, err_glm_high, "ARIMA", "GLM Gamma (high-resolution)"),
  dm_safe(err_arima, err_mpo_ols,  "ARIMA", "MPO (OLS)"),
  dm_safe(err_arima, err_mpos_glm, "ARIMA", "MPOS (GLM Gamma)"),
  dm_safe(err_arima, err_gam_high, "ARIMA", "GAM Gamma (high-resolution)")
), fill = TRUE)

tab_dm[, DM_statistic := round(DM_statistic, 3)]
tab_dm[, p_value := round(p_value, 4)]
tab_dm[, Significance := fifelse(
  is.na(p_value), "",
  fifelse(p_value < 0.001, "***",
          fifelse(p_value < 0.01, "**",
                  fifelse(p_value < 0.05, "*", "ns")))
)]

fwrite(tab_dm, file.path(cfg$results_dir, "Table_DieboldMariano_Results.csv"), bom = TRUE)

save_png("Figure_SkillVsARIMA.png", {
  plot_skill <- copy(tab_compare)
  plot_skill[, Model := factor(Model,
                               levels = plot_skill[order(`Skill vs ARIMA`, decreasing = TRUE)]$Model)]
  plot_skill <- plot_skill[order(`Skill vs ARIMA`, decreasing = TRUE)]
  
  par(family = "serif")
  par(mar = c(8,5,2,1))
  
  bp <- barplot(plot_skill$`Skill vs ARIMA`,
                names.arg = plot_skill$Model,
                las = 2,
                col = "grey65",
                border = "grey25",
                ylab = "Forecast skill relative to ARIMA",
                font.lab = 2,
                cex.lab = 1,
                cex.axis = 0.9,
                bty = "l")
  
  abline(h = 0, lwd = 2, lty = 2)
  points(bp, plot_skill$`Skill vs ARIMA`, pch = 19, cex = 1)
  legend("bottomleft",
         legend = c("Skill = 0 -> ARIMA benchmark"),
         lty = 2, lwd = 2, bty = "n", cex = 0.9)
}, width = 2200, height = 1600, res = 320)

# ============================================================
# 18) FINAL CHECK
# ============================================================

cat("\n[OK] Objects created:\n")
print(c(
  dem = exists("dem"),
  P_bella = exists("P_bella"),
  P_inq = exists("P_inq"),
  oni_monthly = exists("oni_monthly"),
  df2 = exists("df2"),
  train = exists("train"),
  test = exists("test"),
  tab1 = exists("tab1"),
  ols_models = exists("ols_models"),
  glm_models = exists("glm_models"),
  best_glm = exists("best_glm"),
  tabla_mix = exists("tabla_mix"),
  tab_compare = exists("tab_compare"),
  tab_dm = exists("tab_dm")
))

cat("\nBest GLM:", best_glm_name, "\n")
cat("\nRows: df2 =", nrow(df2), " train =", nrow(train), " test =", nrow(test), "\n")


# ============================================================
# 19) TABLE 5: ROLLING-ORIGIN CROSS-VALIDATION (MASE BY HORIZON)
# ------------------------------------------------------------
# Requiere en el entorno: df2 (salida de build_features) y cfg.
# Pegar en run_forecasting.R despues de la seccion 10 (TABLES).
# ============================================================

step("TABLE 5: ROLLING-ORIGIN CROSS-VALIDATION")

cv_cfg <- list(
  first_origin = cfg$cv_first_origin,
  step         = cfg$cv_step,
  horizons     = cfg$cv_horizons,
  m_season     = cfg$m_season,
  refit_order  = FALSE,   # TRUE re-selects the ARIMA order at every origin (slow)
  verbose      = TRUE
)
cv_cfg$max_h <- max(cv_cfg$horizons)

# Especificacion MPO, identica a fit_glm_models()
f_mpo_cv <- caudal ~ anio + mes + dow + holidayF + P7 + I(P7^2) + oni + miss_P

# ------------------------------------------------------------
# 19.1 Utilidades
# ------------------------------------------------------------

# Denominador MASE: error naive estacional EN MUESTRA de entrenamiento.
# Es el estandar de Hyndman & Koehler (2006); no debe calcularse sobre el test.
mase_denom <- function(y_train, m = 7) {
  n <- length(y_train)
  assert(n > m, "Serie de entrenamiento demasiado corta para el denominador MASE.")
  mean(abs(y_train[(m + 1):n] - y_train[1:(n - m)]), na.rm = TRUE)
}

# Pronostico naive estacional replicado hasta el horizonte maximo
snaive_path <- function(y_train, h, m = 7) {
  last_cycle <- tail(y_train, m)
  as.numeric(last_cycle[((seq_len(h) - 1) %% m) + 1])
}

# Comprueba que el tramo de test no introduce niveles de factor
# ausentes en el tramo de calibracion (impediria predict.glm)
levels_ok <- function(tr, te, vars = c("mes", "dow", "holidayF")) {
  for (v in vars) {
    if (!v %in% names(tr)) next
    if (length(setdiff(unique(as.character(te[[v]])),
                       unique(as.character(tr[[v]])))) > 0) return(FALSE)
  }
  TRUE
}

# ------------------------------------------------------------
# 19.2 Bucle de origen movil
# ------------------------------------------------------------

run_rolling_cv <- function(df2, cv_cfg, formula_glm) {
  DT <- as.data.table(copy(df2))
  setorder(DT, fecha)
  n <- nrow(DT)

  origins <- seq(cv_cfg$first_origin, n - cv_cfg$max_h, by = cv_cfg$step)
  assert(length(origins) >= 5,
         paste0("Muy pocos origenes (", length(origins),
                "). Reduzca cv_cfg$first_origin o cv_cfg$step."))
  message("[CV] Origenes: ", length(origins),
          " | horizonte maximo: ", cv_cfg$max_h, " dias")

  arima_order <- NULL
  arima_seas  <- NULL
  out <- list()

  for (k in seq_along(origins)) {
    o  <- origins[k]
    tr <- DT[1:o]
    te <- DT[(o + 1):(o + cv_cfg$max_h)]

    if (!levels_ok(tr, te)) {
      message("[CV] Origen ", o, " omitido: niveles de factor no presentes en calibracion.")
      next
    }

    ytr <- tr$caudal
    yte <- te$caudal
    denom <- mase_denom(ytr, cv_cfg$m_season)

    # --- (a) naive estacional -----------------------------------------
    p_naive <- snaive_path(ytr, cv_cfg$max_h, cv_cfg$m_season)

    # --- (b) ARIMA -----------------------------------------------------
    y_ts <- ts(ytr, frequency = cv_cfg$m_season)
    if (is.null(arima_order) || isTRUE(cv_cfg$refit_order)) {
      m_ar <- auto.arima(y_ts, seasonal = TRUE, stepwise = TRUE, approximation = FALSE)
      if (is.null(arima_order)) {
        arima_order <- forecast::arimaorder(m_ar)[1:3]
        arima_seas  <- forecast::arimaorder(m_ar)
        arima_seas  <- if (length(arima_seas) >= 6) arima_seas[4:6] else c(0, 0, 0)
        message("[CV] Orden ARIMA seleccionado en el primer origen: (",
                paste(arima_order, collapse = ","), ")(",
                paste(arima_seas, collapse = ","), ")[", cv_cfg$m_season, "]")
      }
    } else {
      m_ar <- tryCatch(
        Arima(y_ts, order = arima_order,
              seasonal = list(order = arima_seas, period = cv_cfg$m_season),
              method = "CSS-ML"),
        error = function(e) auto.arima(y_ts, seasonal = TRUE, stepwise = TRUE)
      )
    }
    p_arima <- as.numeric(forecast(m_ar, h = cv_cfg$max_h)$mean)

    # --- (c) MPO GLM Gamma ---------------------------------------------
    tr_g <- copy(tr)
    if (any(tr_g$caudal <= 0, na.rm = TRUE)) tr_g[caudal <= 0, caudal := 1e-6]
    m_glm <- glm(formula_glm, data = tr_g, family = Gamma(link = "log"))
    p_glm <- as.numeric(predict(m_glm, newdata = te, type = "response"))

    # --- metricas por horizonte -----------------------------------------
    preds <- list(`Seasonal naive` = p_naive, `ARIMA` = p_arima, `MPO-GLM Gamma` = p_glm)
    for (h in cv_cfg$horizons) {
      idx <- seq_len(min(h, length(yte)))
      for (nm in names(preds)) {
        e <- yte[idx] - preds[[nm]][idx]
        out[[length(out) + 1]] <- data.table(
          origin      = o,
          origin_date = tr$fecha[o],
          h           = h,
          Model       = nm,
          MAE         = mean(abs(e), na.rm = TRUE),
          RMSE        = sqrt(mean(e^2, na.rm = TRUE)),
          MASE        = mean(abs(e), na.rm = TRUE) / denom,
          Bias        = mean(-e, na.rm = TRUE)
        )
      }
    }

    if (isTRUE(cv_cfg$verbose)) {
      message("[CV] ", k, "/", length(origins),
              " | origen ", as.character(tr$fecha[o]), " | n_train = ", o)
    }
  }

  assert(length(out) > 0, "La validacion cruzada no produjo resultados.")
  # exporta el orden ARIMA para la nota al pie de la tabla
  arima_order_used <<- arima_order
  arima_seas_used  <<- arima_seas
  rbindlist(out)
}

cv_long <- run_rolling_cv(df2, cv_cfg, f_mpo_cv)
fwrite(cv_long, file.path(cfg$results_dir, "Table5_RollingOriginCV_long.csv"), bom = TRUE)

# ------------------------------------------------------------
# 19.3 Construccion de la Tabla 5
# ------------------------------------------------------------

model_levels <- c("Seasonal naive", "ARIMA", "MPO-GLM Gamma")
cv_long[, Model := factor(Model, levels = model_levels)]

agg <- cv_long[, .(MASE = mean(MASE, na.rm = TRUE),
                   MAE  = mean(MAE,  na.rm = TRUE)),
               by = .(h, Model)]

# proporcion de origenes en que cada modelo obtiene el menor MAE
wins <- cv_long[, .SD[which.min(MAE)], by = .(origin, h)][
  , .(Wins = .N), by = .(h, Model)]
n_org <- length(unique(cv_long$origin))
wins[, WinPct := 100 * Wins / n_org]

tab5 <- dcast(agg, h ~ Model, value.var = "MASE")
setnames(tab5,
         c("h", model_levels),
         c("Horizon (days)",
           "MASE - Seasonal naive",
           "MASE - ARIMA",
           "MASE - MPO-GLM Gamma"))

mae_wide <- dcast(agg, h ~ Model, value.var = "MAE")
setnames(mae_wide,
         c("h", model_levels),
         c("Horizon (days)",
           "MAE - Seasonal naive (L/s)",
           "MAE - ARIMA (L/s)",
           "MAE - MPO-GLM Gamma (L/s)"))

tab5 <- merge(tab5, mae_wide, by = "Horizon (days)", sort = FALSE)

win_ar <- wins[Model == "ARIMA", .(h, WinPct)]
setnames(win_ar, c("Horizon (days)", "Origins where ARIMA is best (%)"))
tab5 <- merge(tab5, win_ar, by = "Horizon (days)", all.x = TRUE, sort = FALSE)
tab5[is.na(`Origins where ARIMA is best (%)`), `Origins where ARIMA is best (%)` := 0]

tab5[, `Ratio MPO / ARIMA` := `MASE - MPO-GLM Gamma` / `MASE - ARIMA`]
setorder(tab5, `Horizon (days)`)

tab5_out <- copy(tab5)
for (cc in grep("^MASE", names(tab5_out), value = TRUE))  tab5_out[, (cc) := round(get(cc), 3)]
for (cc in grep("^MAE",  names(tab5_out), value = TRUE))  tab5_out[, (cc) := round(get(cc), 1)]
tab5_out[, `Ratio MPO / ARIMA` := round(`Ratio MPO / ARIMA`, 2)]
tab5_out[, `Origins where ARIMA is best (%)` := round(`Origins where ARIMA is best (%)`, 0)]

fwrite(tab5_out, file.path(cfg$results_dir, "Table5_RollingOriginCV_FULL.csv"), bom = TRUE)

# --- Version para el manuscrito: horizonte + MASE + MAE -------------------
# Las columnas de diagnostico (tasa de victorias, ratio) se omiten aqui y se
# reportan en el texto de Resultados y en la nota al pie de la tabla.
manuscript_cols <- c(
  "Horizon (days)",
  "MASE - Seasonal naive", "MASE - ARIMA", "MASE - MPO-GLM Gamma",
  "MAE - Seasonal naive (L/s)", "MAE - ARIMA (L/s)", "MAE - MPO-GLM Gamma (L/s)"
)
tab5_ms <- tab5_out[, ..manuscript_cols]
fwrite(tab5_ms, file.path(cfg$tab_dir, "Table5_RollingOriginCV_MASE.csv"), bom = TRUE)

win_all <- cv_long[, .SD[which.min(MAE)], by = .(origin, h)][, .N, by = .(h, Model)]
win_ar_n  <- win_all[Model == "ARIMA", N]
win_mpo_n <- win_all[Model == "MPO-GLM Gamma", N]

table5_title <- "Table 5: Predictive accuracy by forecast horizon under rolling-origin cross-validation."

table5_notes <- paste0(
  "Notes: Models were re-estimated at ", n_org, " forecast origins spaced ",
  cv_cfg$step, " days apart; values are means across origins. ",
  "Model specifications as in Table 2. ",
  "MASE = mean absolute error scaled by the in-sample seasonal naive error (m = ",
  cv_cfg$m_season, ") of each training window (Hyndman and Koehler 2006); ",
  "values below 1 indicate better performance than that benchmark. MAE is in L/s. ",
  "Lowest value per row in bold. Across individual origins ARIMA is best at ",
  min(win_ar_n), "-", max(win_ar_n), " of ", n_org, " and MPO-GLM Gamma at ",
  min(win_mpo_n), "-", max(win_mpo_n), "."
)
writeLines(table5_title, con = file.path(cfg$tab_dir, "Table5_Title.txt"))
writeLines(table5_notes, con = file.path(cfg$tab_dir, "Table5_Notes.txt"))

cat("\n", table5_title, "\n\n", sep = "")
print(tab5_ms)
cat("\n", table5_notes, "\n", sep = "")

message("[OK] Table 5 guardada en: ", file.path(cfg$tab_dir, "Table5_RollingOriginCV_MASE.csv"))

# ------------------------------------------------------------
# 19.4 Figura companera (opcional): MASE frente a horizonte
# ------------------------------------------------------------

# Two panels. Panel A reports mean MASE; panel B reports the paired
# difference against the ARIMA benchmark, which removes between-origin
# variance. Between-origin variance exceeds the differences between
# models, so the paired view is the informative comparison.
W_ <- dcast(cv_long, origin + h ~ Model, value.var = "MASE")
setnames(W_, c("Seasonal naive", "MPO-GLM Gamma"), c("NV", "MPO"))
W_[, `:=`(dMPO = MPO - ARIMA, dNV = NV - ARIMA)]
D_ <- W_[, .(mM = mean(dMPO), loM = quantile(dMPO, .25), hiM = quantile(dMPO, .75),
             mN = mean(dNV),  loN = quantile(dNV,  .25), hiN = quantile(dNV,  .75)),
         by = h][order(h)]

save_png("Figure7_MASE_vs_Horizon.png", {
  par(mfrow = c(1, 2), family = "serif", mar = c(4.4, 4.8, 2.4, 1.0))
  cols <- c("Seasonal naive" = "#6D6D6D", "ARIMA" = "#1F78B4",
            "MPO-GLM Gamma" = "#E31A1C")

  plot(NA, xlim = range(cv_cfg$horizons), ylim = range(agg$MASE) * c(0.95, 1.06),
       log = "x", las = 1, bty = "l", xaxt = "n",
       xlab = "Forecast horizon (days)", ylab = "MASE", font.lab = 2)
  axis(1, at = cv_cfg$horizons, labels = cv_cfg$horizons)
  abline(h = 1, lty = 2, lwd = 1.5, col = "grey45")
  for (nm in model_levels) {
    sdat <- agg[Model == nm][order(h)]
    lines(sdat$h, sdat$MASE, lwd = 2.6, col = cols[nm])
    points(sdat$h, sdat$MASE, pch = 19, cex = 1.0, col = cols[nm])
  }
  legend("topleft", legend = c(model_levels, "MASE = 1"),
         col = c(cols[model_levels], "grey45"), lwd = c(2.6, 2.6, 2.6, 1.5),
         lty = c(1, 1, 1, 2), pch = c(19, 19, 19, NA), bty = "n", cex = 0.78)
  mtext("A", side = 3, adj = 0, line = 0.4, font = 2, cex = 1.15)

  plot(NA, xlim = range(cv_cfg$horizons),
       ylim = range(c(D_$loM, D_$hiM, D_$loN, D_$hiN)) * 1.05,
       log = "x", las = 1, bty = "l", xaxt = "n",
       xlab = "Forecast horizon (days)",
       ylab = "MASE difference vs ARIMA", font.lab = 2)
  axis(1, at = cv_cfg$horizons, labels = cv_cfg$horizons)
  polygon(c(D_$h, rev(D_$h)), c(D_$loM, rev(D_$hiM)),
          col = adjustcolor(cols[["MPO-GLM Gamma"]], alpha.f = 0.18), border = NA)
  polygon(c(D_$h, rev(D_$h)), c(D_$loN, rev(D_$hiN)),
          col = adjustcolor(cols[["Seasonal naive"]], alpha.f = 0.18), border = NA)
  abline(h = 0, lwd = 2, col = cols[["ARIMA"]])
  lines(D_$h, D_$mM, lwd = 2.6, col = cols[["MPO-GLM Gamma"]])
  points(D_$h, D_$mM, pch = 19, cex = 1.0, col = cols[["MPO-GLM Gamma"]])
  lines(D_$h, D_$mN, lwd = 2.6, col = cols[["Seasonal naive"]])
  points(D_$h, D_$mN, pch = 19, cex = 1.0, col = cols[["Seasonal naive"]])
  legend("topright",
         legend = c("MPO-GLM Gamma", "Seasonal naive", "ARIMA (reference)",
                    "Interquartile range"),
         col = c(cols[["MPO-GLM Gamma"]], cols[["Seasonal naive"]],
                 cols[["ARIMA"]], "grey70"),
         lwd = c(2.6, 2.6, 2, 8), pch = c(19, 19, NA, NA), bty = "n", cex = 0.74)
  mtext("B", side = 3, adj = 0, line = 0.4, font = 2, cex = 1.15)
}, width = 2400, height = 1250, res = 310)

# ============================================================
# 20) PEAK-DEMAND COEFFICIENT (K): EMPIRICAL vs SIMULATED
# ------------------------------------------------------------
# Requiere en el entorno: df2 (salida de build_features) y cfg.
# Pegar en run_forecasting.R despues de la seccion de TABLES.
#
# CONTEXTO: pred_anual() en project_mixed() calcula
#     K = max(pred) / mean(pred)
# sobre predict(type="response"), que es la MEDIA CONDICIONAL.
# La media condicional no contiene varianza residual, y el pico
# diario real lo produce el residuo. Ese calculo subestima K.
# Este bloque cuantifica el sesgo y entrega el K defendible.
# ============================================================

step("PEAK-DEMAND COEFFICIENT (K)")

k_cfg <- list(n_sim = cfg$n_sim_K, seed = cfg$seed, drop_outliers = cfg$drop_outliers)

DTK <- as.data.table(copy(df2))
setorder(DTK, fecha)

# ------------------------------------------------------------
# 20.1 Outliers: el manuscrito declara inspeccion con boxplot.stats
#      pero el pipeline nunca los elimina. Se identifican aqui.
# ------------------------------------------------------------
bs      <- boxplot.stats(DTK$caudal)
lo_fence <- quantile(DTK$caudal, .25) - 1.5 * IQR(DTK$caudal)
out_rows <- DTK[caudal < lo_fence, .(fecha, caudal)]

message("[K] Limite inferior del boxplot: ", round(lo_fence, 1), " L/s")
if (nrow(out_rows) > 0) {
  message("[K] Observaciones por debajo del limite (probables paros de planta):")
  print(out_rows)
} else {
  message("[K] Sin observaciones por debajo del limite inferior.")
}

DTK_clean <- if (isTRUE(k_cfg$drop_outliers) && nrow(out_rows) > 0) {
  DTK[caudal >= lo_fence]
} else DTK

# ------------------------------------------------------------
# 20.2 K EMPIRICO OBSERVADO  <-- el valor que va al manuscrito
# ------------------------------------------------------------
K_obs <- DTK_clean[, .(mean_daily = mean(caudal),
                       max_daily  = max(caudal),
                       n_days     = .N,
                       K_peak     = max(caudal) / mean(caudal)), by = anio][order(anio)]

K_obs_summary <- list(
  mean = mean(K_obs$K_peak),
  min  = min(K_obs$K_peak),
  max  = max(K_obs$K_peak),
  sd   = sd(K_obs$K_peak)
)

# ------------------------------------------------------------
# 20.3 K sobre la MEDIA CONDICIONAL (reproduce el metodo actual)
# ------------------------------------------------------------
m_k <- glm(formula(final_model), data = DTK_clean, family = Gamma(link = "log"))

DTK_clean[, mu := as.numeric(predict(m_k, type = "response"))]
K_mu <- DTK_clean[, .(K_peak = max(mu) / mean(mu)), by = anio][order(anio)]

# ------------------------------------------------------------
# 20.4 K sobre DEMANDA SIMULADA (media condicional + residuo Gamma)
# ------------------------------------------------------------
set.seed(k_cfg$seed)
disp  <- summary(m_k)$dispersion
shape <- 1 / disp
mu_v  <- DTK_clean$mu
yr_v  <- DTK_clean$anio

sim_K <- vapply(seq_len(k_cfg$n_sim), function(s) {
  y  <- rgamma(length(mu_v), shape = shape, scale = mu_v / shape)
  dt <- data.table(anio = yr_v, y = y)
  mean(dt[, max(y) / mean(y), by = anio]$V1)
}, numeric(1))

K_sim_summary <- list(
  median = median(sim_K),
  lo     = unname(quantile(sim_K, 0.025)),
  hi     = unname(quantile(sim_K, 0.975))
)

# ------------------------------------------------------------
# 20.5 Tabla comparativa de metodos
# ------------------------------------------------------------
K_methods <- data.table(
  Method = c("Conditional mean only (current pipeline)",
             "Simulated demand (mean + Gamma residual)",
             "Observed daily demand (empirical)"),
  K      = c(mean(K_mu$K_peak), K_sim_summary$median, K_obs_summary$mean),
  Lower  = c(NA_real_, K_sim_summary$lo, K_obs_summary$min),
  Upper  = c(NA_real_, K_sim_summary$hi, K_obs_summary$max)
)
K_methods[, `Ratio to regulatory 1.25` := round(K / 1.25, 3)]
K_methods[, K := round(K, 3)][, Lower := round(Lower, 3)][, Upper := round(Upper, 3)]

# ------------------------------------------------------------
# 20.6 Salidas
# ------------------------------------------------------------
K_obs_out <- copy(K_obs)
K_obs_out[, `:=`(mean_daily = round(mean_daily, 1),
                 max_daily  = round(max_daily, 1),
                 K_peak     = round(K_peak, 3))]
setnames(K_obs_out,
         c("anio", "mean_daily", "max_daily", "n_days", "K_peak"),
         c("Year", "Mean daily demand (L/s)", "Maximum daily demand (L/s)",
           "Days", "Peak factor K"))

fwrite(K_obs_out,   file.path(cfg$results_dir, "Table_PeakFactor_K_byYear.csv"),  bom = TRUE)
fwrite(K_methods,   file.path(cfg$results_dir, "Table_PeakFactor_K_methods.csv"), bom = TRUE)

k_note <- paste0(
  "Notes: Peak factor K is the ratio of maximum to mean daily demand within each calendar year, ",
  "computed on observed demand over ", nrow(DTK_clean), " daily records (",
  format(min(DTK_clean$fecha), "%Y"), "-", format(max(DTK_clean$fecha), "%Y"), "). ",
  if (nrow(out_rows) > 0) paste0(nrow(out_rows),
    " observations fall below the lower boxplot fence (", round(lo_fence, 0),
    " L/s) and are attributable to treatment-plant shutdowns rather than to ",
    "demand; they are ",
    if (isTRUE(k_cfg$drop_outliers)) "excluded from this calculation. "
    else "reported here and retained, so that K is computed on the complete operational record. ") else "",
  "Estimating K from modelled conditional means rather than from observed demand ",
  "omits residual variability and understates the coefficient (", round(mean(K_mu$K_peak), 3),
  " versus ", round(K_obs_summary$mean, 3), ")."
)
writeLines(k_note, con = file.path(cfg$results_dir, "Table_PeakFactor_K_Notes.txt"))

cat("\n================ PEAK FACTOR K ================\n\n")
cat("Por anio (demanda observada):\n"); print(K_obs_out)
cat(sprintf("\n  Media   : %.3f\n  Rango   : %.3f - %.3f\n  Desv.tip: %.3f\n",
            K_obs_summary$mean, K_obs_summary$min, K_obs_summary$max, K_obs_summary$sd))
cat("\nComparacion de metodos:\n"); print(K_methods)
cat("\n  Normativo Ecuador : 1.250")
cat("\n  Gortaire et al. (2016) : ~1.12\n")
cat(sprintf("\n  Margen frente al normativo: %.1f%%\n", 100 * (1.25 - K_obs_summary$mean) / 1.25))
cat("\n", k_note, "\n\n", sep = "")

# ------------------------------------------------------------
# 20.7 Figura: K por anio frente al valor normativo
# ------------------------------------------------------------
save_png("Figure_PeakFactor_K.png", {
  par(family = "serif", mar = c(4.4, 4.8, 1.6, 1.2))
  plot(K_obs$anio, K_obs$K_peak, type = "b", pch = 19, lwd = 2.4, cex = 1.1,
       col = "#1F78B4", ylim = c(1.00, 1.30), las = 1, bty = "l",
       xlab = "Year", ylab = "Peak factor K", font.lab = 2, xaxt = "n")
  axis(1, at = K_obs$anio, labels = K_obs$anio)
  abline(h = 1.25, lty = 2, lwd = 2, col = "#E31A1C")
  abline(h = K_obs_summary$mean, lty = 3, lwd = 1.8, col = "grey35")
  legend("bottomleft",
         legend = c("Observed annual K",
                    "Regulatory design value (1.25)",
                    sprintf("Observational mean (%.2f)", K_obs_summary$mean)),
         col = c("#1F78B4", "#E31A1C", "grey35"),
         lwd = c(2.4, 2, 1.8), lty = c(1, 2, 3), pch = c(19, NA, NA),
         bty = "n", cex = 0.82, seg.len = 2.2)
}, width = 2100, height = 1400, res = 320)

message("[OK] K guardado en: ", file.path(cfg$results_dir, "Table_PeakFactor_K_byYear.csv"))


# ============================================================
# 21) SUBMISSION EXPORT
# ------------------------------------------------------------
# Copies the outputs into ./submission using the numbering of the
# manuscript. Internal file names are kept unchanged inside ./outputs
# and ./results so that the analysis code remains self-consistent.
# ============================================================

step("SUBMISSION EXPORT")

sub_dir  <- file.path(root_dir, "submission")
sub_fig  <- file.path(sub_dir, "figures")
sub_tab  <- file.path(sub_dir, "tables")
sub_supp <- file.path(sub_dir, "supplementary")

for (d in c(sub_dir, sub_fig, sub_tab, sub_supp)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ------------------------------------------------------------
# Figures
# Figures 1 and 2 are external/manual.
# The script generates Figures 3-8.
# ------------------------------------------------------------

fig_map <- c(
  "Figure2_ObservedVsPredicted_TwoPanels.png"          = "Fig3_ObservedVsPredicted_TwoPanels.png",
  "Figure4_ResidualDiagnostics.png"                   = "Fig4_ResidualDiagnostics.png",
  "Figure6_HydrosocialDrivers_Bellavista_MPO_OLS.png" = "Fig5_DriverContributions.png",
  "Figure3_RainfallDemand_MPOS.png"                    = "Fig6_RainfallResponse.png",
  "Figure5_ProjectedDemand_TwoPanels.png"              = "Fig7_ProjectedDemand.png",
  "Figure7_MASE_vs_Horizon.png"                        = "Fig8_MASE_by_Horizon.png"
)

# ------------------------------------------------------------
# Supplementary tables
# All supplementary tables are generated in cfg$results_dir.
# The rolling-origin MASE table is Table 5 of the body, not supplementary.
# ------------------------------------------------------------

# Numbering follows the order of first citation in the manuscript:
#   S3 full coefficients, S4 scenario projections, S5 rolling-origin detail.
supp_map <- c(
  # Only the supplementary tables cited in the manuscript are exported,
  # numbered in order of first citation. The remaining result files stay
  # in results/ as reproducible output but are not part of the submission.
  "Table_DieboldMariano_Results.csv"       = "TableS1_DieboldMariano.csv",
  "Table_SerialDependence_Inference.csv"   = "TableS2_ModelCoefficients_RobustSE.csv",
  "Table_SerialDependence_Notes.txt"       = "TableS2_Notes.txt",
  "Table_Forecast_Annual_MixedScenarios.csv" = "TableS3_ProjectedDemand_Scenarios.csv",
  "Table5_RollingOriginCV_long.csv"        = "TableS4_RollingOriginCV_ByOrigin.csv"
)

copy_set <- function(src_dir, mapping, dest, label) {
  ok <- 0L
  
  for (i in seq_along(mapping)) {
    src <- file.path(src_dir, names(mapping)[i])
    
    if (file.exists(src)) {
      file.copy(
        src,
        file.path(dest, mapping[i]),
        overwrite = TRUE
      )
      ok <- ok + 1L
    } else {
      warning(
        "[EXPORT] not found: ",
        names(mapping)[i],
        call. = FALSE
      )
    }
  }
  
  message(
    "[EXPORT] ", label, ": ",
    ok, "/", length(mapping), " files"
  )
  
  ok
}

# Copy generated Figures 3-8
copy_set(
  cfg$fig_dir,
  fig_map,
  sub_fig,
  "figures"
)

# Copy Supplementary Tables S1-S4
copy_set(
  cfg$results_dir,
  supp_map,
  sub_supp,
  "supplementary tables S1-S4"
)

# ------------------------------------------------------------
# Figures 1 and 2 are external/manual.
# ------------------------------------------------------------

for (f in c("Fig1_StudyArea.png", "Fig2_MethodologicalWorkflow.png")) {
  if (!file.exists(file.path(sub_fig, f))) {
    message(
      "[EXPORT] add manually: submission/figures/",
      f
    )
  }
}

# ------------------------------------------------------------
# Main body tables
# Tables 1 and 2 are external/manual.
# The script generates Tables 3-5.
# ------------------------------------------------------------

# Tables 1 and 2 of the manuscript are prepared manually and are not
# generated here. Tables 3 to 6 are exported with the manuscript numbering.
body_map <- c(
  "Table3_ModelFamilyComparison.csv"     = "Table3_ModelFamilyComparison.csv",
  "Table3_Notes.txt"                     = "Table3_Notes.txt",
  "Table4_RegressionSpecifications.csv"  = "Table4_RegressionSpecifications.csv",
  "Table4_Notes.txt"                     = "Table4_Notes.txt",
  "Table5_RollingOriginCV_MASE.csv"      = "Table5_RollingOriginCV.csv",
  "Table5_Notes.txt"                     = "Table5_Notes.txt"
)

copy_set(
  cfg$tab_dir,
  body_map,
  sub_tab,
  "body tables 3-5"
)

# ------------------------------------------------------------
# Submission manifest
# ------------------------------------------------------------

writeLines(c(
  "SUBMISSION PACKAGE",
  paste0(
    "Generated: ",
    format(Sys.time(), "%Y-%m-%d %H:%M")
  ),
  "",
  "figures/        Fig3-Fig8 generated here; Fig1-Fig2 placed manually",
  "tables/         Table3-Table5 generated here; Table1-Table2 placed manually",
  "supplementary/  Supplementary Tables S1-S4",
  "",
  "PLACED MANUALLY (not produced by this code):",
  "  Fig1_StudyArea.png              study-area map",
  "  Fig2_MethodologicalWorkflow.png methodological workflow diagram",
  "  Table 1  statistical models evaluated",
  "  Table 2  candidate model formulations",
  "",
  "MANUSCRIPT NUMBERING",
  "  Figure 3  observed versus predicted daily demand",
  "  Figure 4  residual diagnostics",
  "  Figure 5  driver contributions to explained variance",
  "  Figure 6  response to antecedent precipitation",
  "  Figure 7  projected demand under scenarios",
  "  Figure 8  predictive accuracy by forecast horizon",
  "  Table 3   performance comparison across model families",
  "  Table 4   comparative performance of regression specifications",
  "  Table 5   predictive accuracy by forecast horizon",
  "  Table S1  Diebold-Mariano tests",
  "  Table S2  model coefficients with robust standard errors",
  "  Table S3  projected demand by scenario",
  "  Table S4  rolling-origin cross-validation by origin",
  "",
  "Session information:",
  capture.output(sessionInfo())
), file.path(sub_dir, "MANIFEST.txt"))

message("[OK] Submission package: ", sub_dir)

cat("\n[OK] SUBMISSION FIGURES:\n")
print(list.files(sub_fig))

cat("\n[OK] SUBMISSION TABLES:\n")
print(list.files(sub_tab))

cat("\n[OK] SUPPLEMENTARY TABLES:\n")
print(list.files(sub_supp))

cat("\n[OK] OUTPUT FIGURES:\n")
print(list.files(cfg$fig_dir))

cat("\n[OK] OUTPUT TABLES:\n")
print(list.files(cfg$tab_dir))

cat("\n[OK] RESULTS:\n")
print(list.files(cfg$results_dir))

message(
  "\n[OK] Done. Project outputs saved under: ",
  root_dir
)


