# ============================================================
# 00_predictive_decomposition_nested_cv.R
#
# Reproducible analysis for:
#   "Wearable sensing reveals the structure of cardiac
#    activation associated with everyday driving"
#
# PURPOSE
#   Fit and evaluate the manuscript's predictive decomposition of
#   instantaneous heart rate separately in two low-motion contexts:
#
#     1. DRIVING
#     2. NONDRIVING_SEDENTARY
#
# SCIENTIFIC MODEL
#   For participant i, time t, day d, and context c:
#
#     HR_raw,itd = HR_base,i + tau_c + phi(X_itd) + epsilon_itd
#
#   Equivalently:
#
#     NHR_itd = HR_raw,itd - HR_base,i
#             = tau_c + phi(X_itd) + epsilon_itd
#
#   HR_base,i is the participant-level mean of the available daily
#   Apple HealthKit resting-heart-rate values. tau_c is the mean
#   empirical elevation above participant baseline within the
#   corresponding behavioral stratum. It is a context-associated
#   offset, not a causal effect.
#
# MODELS
#   The script compares three nested prediction stages:
#
#     baseline_only
#       HR_hat = HR_base,i
#
#     baseline_plus_offset
#       HR_hat = HR_base,i + tau_c
#
#     enet
#       HR_hat = HR_base,i + tau_c + learned covariate modulation
#
# VALIDATION
#   - Outer grouped 5-fold cross-validation estimates held-out
#     performance.
#   - Inner grouped 4-fold cross-validation tunes Elastic Net.
#   - Participant is the grouping unit in both loops; therefore, no
#     participant contributes observations to both training and test
#     data within a fold.
#   - The context offset and all data-dependent preprocessing are
#     estimated from training data only.
#
# PRIMARY METRICS
#   - RMSE (bpm)
#   - MAE (bpm; retained as a supplementary diagnostic)
#   - squared Pearson correlation (r^2)
#   - conventional held-out R^2_OOS = 1 - SSE/SST
#
# INPUT
#   Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
# REQUIRED ACTIVITY LABELS
#   activity3 == "driving"
#   activity3 == "non_driving_sedentary"
#
# OUTPUT
#   Results/nubi_ml/<timestamp>_<resolution>sec_predictive_decomposition_nested_cv/
#
#   Principal downstream files include:
#     predictions_all_models_both_strata.csv
#     compare_metrics_rawhr_by_fold_by_stratum.csv
#     compare_metrics_rawhr_pooled_oof_by_stratum.csv
#     manuscript_table4_metrics.csv
#     feature_importance_terms_<STRATUM>.csv
#     feature_importance_grouped_<STRATUM>.csv
#     feature_importance_grouped_modulators_only_<STRATUM>.csv
#     best_params_by_stratum.csv
#
#   Fold-level metrics support manuscript Figure 5 and Table 4.
#   Modulator-only importance outputs support Figure 6, while the
#   saved full-data ENet refits support downstream Figure 7 analyses.
#
# REPOSITORY SCOPE AND PRIVACY
#   The public repository begins with the curated, analysis-ready
#   minute-level dataset; it does not reconstruct that dataset from
#   raw wearable, smartphone, vehicle, or annotation streams. Direct
#   GPS coordinates are neither required nor used by this analysis.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  
  library(tidymodels)
  library(doParallel)
  
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(rlang)
  library(stringr)
  library(scales)
  library(forcats)
})

options(warn = 1)
tidymodels_prefer()

SCRIPT_VERSION <- "1.0.1"
ANALYSIS_SEED <- 20260309L
set.seed(ANALYSIS_SEED)

# ----------------------------
# User toggles
# ----------------------------
LOCAL_TZ <- "America/Chicago"

USE_COMMON_SUBJECTS <- TRUE
V_OUTER <- 5
V_INNER <- 4

# Driving dynamics windows (minutes)
DYN_WINDOWS_MIN <- c(1, 3, 5)

# Applicability pruning thresholds
KEEP_MAX_NA   <- 0.95
KEEP_MIN_UNIQ <- 2

# If TRUE, include trip_id when it is not extremely high-cardinality
ALLOW_TRIP_ID_IF_REASONABLE <- FALSE
TRIP_ID_MAX_LEVEL_SHARE <- 0.50   # require n_unique_trip_id <= 50% of n_rows

# MASTER strata labels
ACT3_DRIVING <- "driving"
ACT3_ND_SED  <- "non_driving_sedentary"

# ----------------------------
# Robust wd = Scripts/
# ----------------------------
this_script <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)

script_dir <- if (!is.na(this_script) && file.exists(this_script)) {
  dirname(this_script)
} else {
  getwd()
}
script_dir <- normalizePath(script_dir, mustWork = TRUE)

setwd(script_dir)

project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

results_root <- file.path(project_root, "Results")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
results_root <- normalizePath(results_root, mustWork = TRUE)

data_root    <- normalizePath(file.path(project_root, "Data"), mustWork = TRUE)

message("Working directory (Scripts): ", getwd())
message("Project root: ", project_root)

# ----------------------------
# Resolution picker
# ----------------------------
# The public repository currently includes the curated 60-second dataset:
#
#   Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
# The picker is retained so the same script can be reused if 10-sec or
# 30-sec analysis datasets are added later.

pick_resolution <- function(default = 60L) {
  cat(
    "\nChoose dataset resolution:\n",
    "  1) 10 sec  [requires Data/NUBI_Data_10sec_Level_MASTER_CLEAN.csv]\n",
    "  2) 30 sec  [requires Data/NUBI_Data_30sec_Level_MASTER_CLEAN.csv]\n",
    "  3) 60 sec  [included in this repository]\n",
    sep = ""
  )
  
  ans <- trimws(readline(
    sprintf("Enter 10 / 30 / 60 (or 1/2/3). Press Enter for %d sec: ", default)
  ))
  
  if (ans == "") return(as.integer(default))
  if (ans %in% c("1", "10")) return(10L)
  if (ans %in% c("2", "30")) return(30L)
  if (ans %in% c("3", "60")) return(60L)
  
  stop("Invalid entry: ", ans, " (expected 10/30/60 or 1/2/3)")
}

RES_SECONDS <- pick_resolution(default = 60L)

in_path <- file.path(
  data_root,
  sprintf("NUBI_Data_%dsec_Level_MASTER_CLEAN.csv", RES_SECONDS)
)

if (!file.exists(in_path)) {
  stop(
    "Dataset not found: ", in_path, "\n",
    "The public repository currently includes only the 60-second dataset. ",
    "Use 60 sec, or add the corresponding ", RES_SECONDS,
    "-second dataset under Data/."
  )
}

# ----------------------------
# Output folder
# ----------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(
  results_root,
  "nubi_ml",
  paste0(stamp, "_", RES_SECONDS, "sec_predictive_decomposition_nested_cv")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(out_dir)) {
  stop("Could not create output directory: ", out_dir)
}

fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(fig_dir)) {
  stop("Could not create figure directory: ", fig_dir)
}

log_file <- file.path(out_dir, "run_log.txt")

# Logging is intentionally non-fatal. Messages are always printed to the
# console; failure to append to run_log.txt must never terminate the analysis.
log_msg <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste(..., collapse = "")
  )

  message(msg)

  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }

  tryCatch(
    cat(msg, "\n", file = log_file, append = TRUE),
    error = function(e) {
      warning(
        "Could not append to log file '", log_file,
        "': ", conditionMessage(e),
        ". Analysis will continue."
      )
      invisible(NULL)
    }
  )
}

log_msg("Script: 00_predictive_decomposition_nested_cv.R")
log_msg("Resolution: ", RES_SECONDS, " sec")
log_msg("Input: ", in_path)
log_msg("Output dir: ", out_dir)
log_msg("Figure dir: ", fig_dir)
log_msg("LOCAL_TZ: ", LOCAL_TZ)
log_msg("USE_COMMON_SUBJECTS: ", USE_COMMON_SUBJECTS)
log_msg("V_OUTER: ", V_OUTER)
log_msg("V_INNER: ", V_INNER)
log_msg("ENet mode: nested grouped CV; inner folds tune hyperparameters, outer folds estimate performance")

# Record the run configuration before analysis begins.
run_config <- tibble(
  parameter = c(
    "script", "script_version", "seed", "local_timezone", "resolution_seconds",
    "common_subjects_only", "outer_folds", "inner_folds",
    "dynamic_windows_minutes", "max_predictor_missingness",
    "minimum_unique_values", "allow_trip_id", "logging_mode"
  ),
  value = c(
    "00_predictive_decomposition_nested_cv.R", SCRIPT_VERSION,
    as.character(ANALYSIS_SEED), LOCAL_TZ,
    as.character(RES_SECONDS), as.character(USE_COMMON_SUBJECTS),
    as.character(V_OUTER), as.character(V_INNER),
    paste(DYN_WINDOWS_MIN, collapse = ","), as.character(KEEP_MAX_NA),
    as.character(KEEP_MIN_UNIQ), as.character(ALLOW_TRIP_ID_IF_REASONABLE),
    "console_always; file_append_nonfatal"
  )
)
write_csv(run_config, file.path(out_dir, "run_configuration.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))

# ============================================================
# Helpers
# ============================================================
snakeify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_time_local <- function(x, tz = LOCAL_TZ) {
  if (inherits(x, "POSIXt")) return(with_tz(x, tzone = tz))
  if (!is.character(x)) x <- as.character(x)
  
  # Important: treat strings as local clock time
  # If there is a trailing Z, strip it rather than interpreting as UTC.
  x2 <- gsub("Z$", "", x)
  
  tt <- suppressWarnings(lubridate::ymd_hms(x2, tz = tz, quiet = TRUE))
  if (!all(is.na(tt))) return(tt)
  
  suppressWarnings(lubridate::parse_date_time(
    x2,
    orders = c("ymd HMS", "ymd HM", "mdy HMS", "mdy HM", "dmy HMS", "dmy HM"),
    tz = tz
  ))
}

canonicalize_names <- function(dt) {
  stopifnot(is.data.table(dt))
  
  old <- names(dt)
  sn  <- snakeify(old)
  
  canon <- c(
    p_id="p_id", pid="p_id", participant_id="p_id", participantid="p_id",
    time="time", timestamp="time", datetime="time", date_time="time",
    
    raw_hr="raw_hr",
    nhr="nhr",
    
    activity="activity",
    activity3="activity3",
    expert_pa="expert_pa",
    
    data_source="data_source", datasource="data_source", source="data_source",
    
    trip_time="trip_time",
    triptime="trip_time",
    trip_id="trip_id",
    trips="trips",
    
    day_num="day_num",
    daynum="day_num",
    days="days",
    day_period="day_period",
    day_type="day_type",
    
    weather_info="weather_info",
    weather="weather_info",
    in_radius="in_radius",
    
    bl_hr="bl_hr",
    hr_bl="bl_hr",
    hrbl="bl_hr",
    baseline_hr="bl_hr",
    baseline="bl_hr",
    hr_baseline="bl_hr",
    
    nr_hr_sd="nr_hr_sd",
    nr_hr_2sd="nr_hr_2sd",
    
    speed="speed",
    atp="atp",
    jf="jf",
    ff="ff",
    ff_speed="ff_speed",
    rtp="rtp",
    energy_acc="energy_acc",
    energy_rot="energy_rot",
    distance="distance",
    trip_distance="trip_distance",
    trip_duration="trip_duration",
    
    live_lat="live_lat",
    live_long="live_long",
    live_lon="live_long",
    home_lat="home_lat",
    home_long="home_long",
    home_lon="home_long",
    first_point_lat="first_point_lat",
    first_point_lon="first_point_lon",
    first_point_long="first_point_lon",
    
    gender="gender",
    age="age",
    trait_anxiety="trait_anxiety",
    morning_anxiety="morning_anxiety",
    
    openness="openness",
    neuroticism="neuroticism",
    conscientiousness="conscientiousness",
    agreeableness="agreeableness",
    extraversion="extraversion",
    
    md="md", pd="pd", td="td", p="p", e="e", f="f"
  )
  
  new <- sn
  hit <- sn %in% names(canon)
  new[hit] <- unname(canon[sn[hit]])
  
  if (!identical(old, new)) {
    new <- make.unique(new, sep = "_")
    setnames(dt, old, new)
  }
  dt
}

safe_save_pdf <- function(plot_obj, path, w = 7, h = 5) {
  ok <- tryCatch({
    ggsave(filename = path, plot = plot_obj, width = w, height = h,
           device = grDevices::cairo_pdf)
    TRUE
  }, error = function(e) FALSE)
  
  if (!ok) {
    ggsave(filename = path, plot = plot_obj, width = w, height = h,
           device = "pdf", useDingbats = FALSE)
  }
}

safe_save_png <- function(plot_obj, path, w = 7, h = 5, dpi = 300) {
  ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}

make_day_key <- function(dt) {
  if ("day_num" %in% names(dt)) {
    x <- dt[["day_num"]]
    if (is.numeric(x) || is.integer(x)) return(as.integer(x))
    xs <- as.character(x)
    dig <- suppressWarnings(as.integer(stringr::str_extract(xs, "\\d+")))
    if (!all(is.na(dig))) return(dig)
  }
  
  if ("days" %in% names(dt)) {
    x <- dt[["days"]]
    if (is.numeric(x) || is.integer(x)) return(as.integer(x))
    xs <- as.character(x)
    dig <- suppressWarnings(as.integer(stringr::str_extract(xs, "\\d+")))
    if (!all(is.na(dig))) return(dig)
  }
  
  if ("dt_time" %in% names(dt) && inherits(dt[["dt_time"]], "POSIXt")) {
    d <- as.Date(dt[["dt_time"]], tz = LOCAL_TZ)
    u <- sort(unique(d))
    return(as.integer(match(d, u)))
  }
  
  NULL
}

compute_metrics_vec <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]

  if (length(y) < 5) {
    return(tibble(
      rmse = NA_real_,
      mae = NA_real_,
      rsq_cor = NA_real_,
      rsq_oos = NA_real_,
      n = length(y)
    ))
  }

  resid <- y - p
  rmse_v <- sqrt(mean(resid^2))
  mae_v  <- mean(abs(resid))

  # Squared-correlation R^2. This is the value reported by
  # yardstick::rsq() and used in the original manuscript.
  rsq_cor_v <- if (sd(y) == 0 || sd(p) == 0) {
    NA_real_
  } else {
    cor(y, p)^2
  }

  # Standard out-of-sample R^2:
  #
  #   1 - sum((y - yhat)^2) / sum((y - mean(y_test))^2)
  #
  # For fold-level summaries, mean(y_test) is the mean of the
  # held-out observations in that fold. For overall out-of-fold
  # summaries, mean(y_test) is the mean across all held-out
  # predictions being summarized.
  sse <- sum(resid^2)
  sst <- sum((y - mean(y))^2)
  rsq_oos_v <- if (sst <= .Machine$double.eps) NA_real_ else 1 - (sse / sst)

  tibble(
    rmse = rmse_v,
    mae = mae_v,
    rsq_cor = rsq_cor_v,
    rsq_oos = rsq_oos_v,
    n = length(y)
  )
}

metrics_to_long <- function(mm) {
  tibble(
    .metric = c("rmse", "mae", "rsq_cor", "rsq_oos"),
    .estimate = c(mm$rmse, mm$mae, mm$rsq_cor, mm$rsq_oos),
    n = mm$n
  )
}


predictor_applicability <- function(df, pred_cols) {
  tibble(var = pred_cols) %>%
    mutate(
      na_rate = map_dbl(var, ~{
        x <- df[[.x]]
        if (is.character(x) || is.factor(x) || is.logical(x)) {
          mean(is.na(x), na.rm = FALSE)
        } else {
          xx <- suppressWarnings(as.numeric(x))
          mean(is.na(xx) | !is.finite(xx), na.rm = FALSE)
        }
      }),
      uniq_n = map_int(var, ~{
        x <- df[[.x]]
        if (is.character(x) || is.factor(x) || is.logical(x)) return(length(unique(x)))
        xx <- suppressWarnings(as.numeric(x))
        length(unique(xx[is.finite(xx)]))
      }),
      class = map_chr(var, ~paste(class(df[[.x]]), collapse = "|"))
    ) %>%
    arrange(desc(na_rate), uniq_n)
}

# ============================================================
# HARD leakage guards (allow ONLY bl_hr_person)
# ============================================================
FORBIDDEN_PRED_REGEX <- paste0(
  "(",
  "^raw_hr($|_)", "|",
  "^nhr($|_)", "|",
  "^hr_bl($|_)", "|",
  "^bl_hr_num($|_)", "|",
  "^bl_hr($|_(?!person$))", "|",
  "^nr_hr_sd($|_)", "|",
  "^nr_hr_2sd($|_)", "|",
  "(_obs$|_hat$|\\.pred$|^\\.pred$|^pred$|^prediction$)",
  ")"
)

assert_no_forbidden_predictors <- function(predictor_names, context = "") {
  bad <- predictor_names[grepl(FORBIDDEN_PRED_REGEX, predictor_names, perl = TRUE)]
  if (length(bad) > 0) {
    stop(
      "LEAKAGE GUARD TRIGGERED", if (nzchar(context)) paste0(" [", context, "]") else "",
      ": forbidden predictors present: ",
      paste(unique(bad), collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_no_forbidden_in_baked <- function(prepped_recipe, te_df, outcome_name = "raw_hr", context = "") {
  baked <- bake(prepped_recipe, new_data = te_df)
  cols_pred <- setdiff(colnames(baked), outcome_name)
  bad <- cols_pred[grepl(FORBIDDEN_PRED_REGEX, cols_pred, perl = TRUE)]
  if (length(bad) > 0) {
    stop(
      "LEAKAGE GUARD TRIGGERED", if (nzchar(context)) paste0(" [", context, "]") else "",
      ": forbidden columns present in baked predictors: ",
      paste(unique(bad), collapse = ", ")
    )
  }
  invisible(TRUE)
}

# ============================================================
# Driving-specific short-term dynamics
# ============================================================
dyn_base_vars <- c(
  "speed", "ff", "ff_speed", "atp", "rtp", "jf",
  "energy_acc", "energy_rot"
)

add_dynamics <- function(dt,
                         id_col = "p_id",
                         time_col = "dt_time",
                         vars,
                         res_seconds,
                         windows_min = c(1, 3, 5)) {
  stopifnot(is.data.table(dt))
  stopifnot(id_col %in% names(dt), time_col %in% names(dt))
  stopifnot(all(vars %in% names(dt)))
  
  # Match the manuscript analysis: dynamics are computed after chronological
  # ordering within participant. Trip ID is not used as a predictor by default.
  group_cols <- id_col
  setorderv(dt, c(id_col, time_col))

  for (v in vars) {
    dt[, (v) := suppressWarnings(as.numeric(get(v)))]
    
    dt[, paste0(v, "_lag1")  := shift(get(v), 1L, type = "lag"), by = group_cols]
    dt[, paste0(v, "_diff1") := get(v) - get(paste0(v, "_lag1")), by = group_cols]
    
    for (wm in windows_min) {
      k <- max(2L, as.integer(round((wm * 60) / res_seconds)))
      tag <- paste0("_", wm, "m")
      
      rm_name <- paste0(v, "_rm", tag)
      rs_name <- paste0(v, "_rs", tag)
      sl_name <- paste0(v, "_slope", tag)
      
      dt[, (rm_name) := frollmean(get(v), n = k, align = "right", fill = NA_real_), by = group_cols]
      dt[, (rs_name) := {
        m1 <- frollmean(get(v), n = k, align = "right", fill = NA_real_)
        m2 <- frollmean(get(v)^2, n = k, align = "right", fill = NA_real_)
        s2 <- m2 - m1^2
        sqrt(pmax(s2, 0))
      }, by = group_cols]
      dt[, (sl_name) := (get(v) - shift(get(v), k - 1L, type = "lag")) /
           ((k - 1L) * res_seconds),
         by = group_cols]
    }
  }
  
  dt
}

# ============================================================
# Feature importance (ENet): |standardized coefficients|
#
# Complete tables retain every fitted predictor, including the
# participant baseline. Separate modulator-only outputs exclude
# bl_hr_person because manuscript Figure 6 concerns covariates
# beyond participant baseline and context-specific offset.
# ============================================================
extract_enet_importance <- function(fitted_wf, penalty_value, predictors_final) {
  eng <- workflows::extract_fit_engine(fitted_wf)
  
  cc <- tryCatch({
    as.matrix(stats::coef(eng, s = penalty_value))
  }, error = function(e) {
    as.matrix(stats::coef(eng))
  })
  
  coef_tbl <- tibble(
    term = rownames(cc),
    estimate = as.numeric(cc[, 1])
  ) %>%
    filter(term != "(Intercept)", is.finite(estimate)) %>%
    mutate(importance = abs(estimate)) %>%
    arrange(desc(importance))
  
  predictors_final <- unique(predictors_final)
  
  coef_tbl <- coef_tbl %>%
    mutate(group = purrr::map_chr(term, function(tt) {
      if (tt %in% predictors_final) return(tt)
      for (p in predictors_final) {
        if (startsWith(tt, paste0(p, "_"))) return(p)
      }
      tt
    }))
  
  group_tbl <- coef_tbl %>%
    group_by(group) %>%
    summarise(
      importance = sum(importance, na.rm = TRUE),
      n_terms = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(desc(importance))
  
  list(term_tbl = coef_tbl, group_tbl = group_tbl)
}

plot_importance_bar <- function(tbl, name_col, value_col = "importance",
                                title = "", top_n = 25) {
  stopifnot(name_col %in% names(tbl), value_col %in% names(tbl))
  
  dd <- tbl %>%
    arrange(desc(.data[[value_col]])) %>%
    slice_head(n = top_n) %>%
    mutate(name = .data[[name_col]]) %>%
    mutate(name = forcats::fct_reorder(name, .data[[value_col]]))
  
  ggplot(dd, aes(x = name, y = .data[[value_col]])) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = "|standardized coefficient| (sum if grouped)", title = title) +
    theme_minimal(base_size = 11)
}

# ============================================================
# Read + standardize schema
# ============================================================
log_msg("Reading FULL data ...")
dt_full <- fread(in_path, showProgress = TRUE)
dt_full <- canonicalize_names(dt_full)

need <- c("p_id", "time", "raw_hr", "activity3", "bl_hr")
miss <- setdiff(need, names(dt_full))
if (length(miss) > 0) {
  stop("Missing required columns after standardization: ", paste(miss, collapse = ", "))
}

dt_full[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]
dt_full <- dt_full[!is.na(dt_time)]
log_msg("Rows after time parse: ", nrow(dt_full), " | subjects: ", uniqueN(dt_full$p_id))

dt_full[, activity3_norm := trimws(tolower(as.character(activity3)))]
act3_levels <- sort(unique(dt_full$activity3_norm))
write_csv(tibble(activity3_level = act3_levels), file.path(out_dir, "activity3_levels.csv"))
log_msg("Wrote activity3_levels.csv (n=", length(act3_levels), "). Levels: ", paste(act3_levels, collapse = ", "))

act3_counts <- dt_full[, .(n_rows = .N, n_subj = uniqueN(p_id)), by = .(activity3_norm)][order(-n_rows)]
write_csv(as.data.frame(act3_counts), file.path(out_dir, "activity3_counts.csv"))
log_msg("Wrote activity3_counts.csv")

# ============================================================
# Participant-level baseline from available participant-day baseline values
# ============================================================
dt_full[, bl_hr_num := suppressWarnings(as.numeric(bl_hr))]
day_key <- make_day_key(dt_full)
if (is.null(day_key)) stop("Could not construct day key. Need day_num/days or dt_time.")
dt_full[, day_key := day_key]

# The curated minute-level file repeats the same participant-day baseline
# across rows. Taking the within-day median defensively recovers one value per
# participant-day while tolerating accidental duplicate-row inconsistencies.
baseline_by_subj_day <- dt_full[
  is.finite(bl_hr_num) & !is.na(day_key),
  .(
    bl_hr_day = median(bl_hr_num, na.rm = TRUE),
    n_rows_used = .N,
    n_unique_bl_hr = uniqueN(bl_hr_num)
  ),
  by = .(p_id, day_key)
]
write_csv(as.data.frame(baseline_by_subj_day), file.path(out_dir, "baseline_by_subject_day.csv"))
log_msg("Wrote baseline_by_subject_day.csv: rows=", nrow(baseline_by_subj_day))
if (any(baseline_by_subj_day$n_unique_bl_hr > 1L)) {
  warning("Some participant-days contain more than one finite baseline value; the within-day median was used.")
}
if (nrow(baseline_by_subj_day) == 0) stop("baseline_by_subject_day has 0 rows. Check bl_hr values.")

baseline_by_subj <- baseline_by_subj_day[
  ,
  .(
    bl_hr_person = mean(bl_hr_day, na.rm = TRUE),
    n_days_with_bl = .N,
    bl_hr_day_sd = sd(bl_hr_day, na.rm = TRUE)
  ),
  by = .(p_id)
]
write_csv(as.data.frame(baseline_by_subj), file.path(out_dir, "baseline_by_subject.csv"))
log_msg("Wrote baseline_by_subject.csv: n_subjects=", nrow(baseline_by_subj))

dt_full <- merge(
  dt_full,
  baseline_by_subj[, .(p_id, bl_hr_person, n_days_with_bl)],
  by = "p_id",
  all.x = TRUE
)
dt_full <- dt_full[is.finite(bl_hr_person)]
log_msg("After dropping missing-baseline subjects: rows=", nrow(dt_full), " | subjects=", uniqueN(dt_full$p_id))

# ============================================================
# Split into strata using activity3
# ============================================================
dt_drive <- dt_full[activity3_norm == ACT3_DRIVING]
dt_nond  <- dt_full[activity3_norm == ACT3_ND_SED]

log_msg("DRIVING rows: ", nrow(dt_drive), " | subjects: ", uniqueN(dt_drive$p_id))
log_msg("NONDRIVING_SEDENTARY rows: ", nrow(dt_nond), " | subjects: ", uniqueN(dt_nond$p_id))

if (nrow(dt_drive) == 0) stop("No DRIVING rows. Check activity3 levels.")
if (nrow(dt_nond) == 0) stop("No NONDRIVING_SEDENTARY rows. Check activity3 levels.")

if (USE_COMMON_SUBJECTS) {
  common_subj <- intersect(unique(dt_drive$p_id), unique(dt_nond$p_id))
  dt_drive <- dt_drive[p_id %in% common_subj]
  dt_nond  <- dt_nond[p_id %in% common_subj]
  
  log_msg("COMMON subjects enforced: ", length(common_subj))
  log_msg("DRIVING rows (common): ", nrow(dt_drive), " | subjects: ", uniqueN(dt_drive$p_id))
  log_msg("NONDRIVING rows (common): ", nrow(dt_nond), " | subjects: ", uniqueN(dt_nond$p_id))
  
  if (uniqueN(dt_drive$p_id) < V_OUTER || uniqueN(dt_nond$p_id) < V_OUTER) {
    stop("After enforcing common subjects, too few subjects for CV. Consider USE_COMMON_SUBJECTS <- FALSE.")
  }
}

# ============================================================
# Add dynamics to DRIVING ONLY
# ============================================================
add_dynamics_if_possible <- function(dt, label) {
  have <- intersect(dyn_base_vars, names(dt))
  if (length(have) == 0) {
    log_msg(label, ": No base vars for dynamics present; skipping dynamics.")
    return(dt)
  }
  
  log_msg(label, ": Adding dynamics for: ", paste(have, collapse = ", "))
  dt <- add_dynamics(
    dt,
    id_col = "p_id",
    time_col = "dt_time",
    vars = have,
    res_seconds = RES_SECONDS,
    windows_min = DYN_WINDOWS_MIN
  )
  log_msg(label, ": Dynamics added.")
  dt
}

dt_drive <- add_dynamics_if_possible(dt_drive, "DRIVING")
log_msg("NONDRIVING_SEDENTARY: dynamics intentionally skipped.")

# ============================================================
# Parallel
# ============================================================
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)
log_msg("Parallel workers: ", n_cores)

# ============================================================
# Predictor selection per stratum
# ============================================================
base_predictors_common <- c(
  "day_num", "days", "day_period", "day_type",
  "trip_time", "in_radius",
  "gender", "age", "trait_anxiety", "morning_anxiety",
  "openness", "neuroticism", "conscientiousness",
  "agreeableness", "extraversion",
  "bl_hr_person"
)

driving_only_core <- c(
  "weather_info",
  "speed", "atp", "jf", "ff", "ff_speed", "rtp",
  "energy_acc", "energy_rot",
  "trip_distance", "trip_duration",
  "md", "pd", "td", "p", "e", "f"
)

select_predictors_for_stratum <- function(dt_stratum, stratum_name, out_dir) {
  cand <- base_predictors_common
  
  if (stratum_name == "DRIVING") {
    cand <- unique(c(cand, driving_only_core))
    
    if (ALLOW_TRIP_ID_IF_REASONABLE && "trip_id" %in% names(dt_stratum)) {
      trip_n_unique <- uniqueN(as.character(dt_stratum$trip_id))
      trip_share <- trip_n_unique / max(1, nrow(dt_stratum))
      if (is.finite(trip_share) && trip_share <= TRIP_ID_MAX_LEVEL_SHARE) {
        cand <- unique(c(cand, "trip_id"))
        log_msg("DRIVING: trip_id allowed as predictor (unique share=", signif(trip_share, 4), ")")
      } else {
        log_msg("DRIVING: trip_id excluded for high cardinality (unique share=", signif(trip_share, 4), ")")
      }
    }
    
    have_dyn <- intersect(dyn_base_vars, names(dt_stratum))
    dyn_cols <- if (length(have_dyn) > 0) {
      grep(
        pattern = paste0("^(", paste(have_dyn, collapse = "|"), ")_"),
        x = names(dt_stratum),
        value = TRUE
      )
    } else {
      character(0)
    }
    cand <- unique(c(cand, dyn_cols))
  }
  
  drop_exact <- c(
    "nr_hr_sd", "nr_hr_2sd",
    "dt_time", "time",
    "data_source", "activity", "activity3", "activity3_norm", "expert_pa",
    "day_key", "n_days_with_bl",
    "bl_hr", "bl_hr_num", "nhr",
    "distance",
    "live_lat", "live_long", "home_lat", "home_long", "first_point_lat", "first_point_lon"
  )
  
  predictors0 <- intersect(cand, names(dt_stratum))
  predictors0 <- setdiff(predictors0, drop_exact)
  predictors0 <- predictors0[!grepl(FORBIDDEN_PRED_REGEX, predictors0, perl = TRUE)]
  
  if (stratum_name != "DRIVING") {
    driving_regex <- "^(md|pd|td|p$|e$|f$|speed|atp|jf|ff|ff_speed|rtp|energy_acc|energy_rot|trip_distance|trip_duration|trip_id|weather_info)($|_)"
    predictors0 <- predictors0[!grepl(driving_regex, predictors0)]
    predictors0 <- setdiff(predictors0, "weather_info")
  }
  
  df_tmp <- dt_stratum[, c("p_id", "raw_hr", predictors0), with = FALSE] |> as.data.frame()
  
  appl <- predictor_applicability(df_tmp, predictors0) %>%
    mutate(stratum = stratum_name, .before = 1)
  
  write_csv(appl, file.path(out_dir, paste0("predictor_applicability_raw_", stratum_name, ".csv")))
  
  good_vars <- appl %>%
    filter(na_rate <= KEEP_MAX_NA, uniq_n >= KEEP_MIN_UNIQ) %>%
    pull(var)
  
  predictors <- good_vars
  assert_no_forbidden_predictors(predictors, context = paste0("predictor_selection_", stratum_name))
  
  out_name <- paste0("predictor_list_", stratum_name, ".csv")
  write_csv(tibble(stratum = stratum_name, predictor = predictors), file.path(out_dir, out_name))
  log_msg(stratum_name, ": FINAL predictor count = ", length(predictors), " (written to ", out_name, ")")
  
  predictors
}

log_msg("Selecting predictors for DRIVING ...")
preds_drive <- select_predictors_for_stratum(dt_drive, "DRIVING", out_dir)

log_msg("Selecting predictors for NONDRIVING_SEDENTARY ...")
preds_nond <- select_predictors_for_stratum(dt_nond, "NONDRIVING_SEDENTARY", out_dir)

write_csv(
  bind_rows(
    tibble(stratum = "DRIVING", predictor = preds_drive),
    tibble(stratum = "NONDRIVING_SEDENTARY", predictor = preds_nond)
  ),
  file.path(out_dir, "predictor_list_BOTH.csv")
)

# ============================================================
# Core evaluation for one stratum
# ============================================================
run_stratum <- function(dt_stratum, stratum_name, predictors_final) {
  log_msg("---- Stratum START: ", stratum_name, " ----")
  
  id_col <- "p_id"
  target_raw <- "raw_hr"
  
  assert_no_forbidden_predictors(predictors_final, context = paste0("run_stratum_", stratum_name))
  
  id_metadata <- intersect(c("dt_time", "day_key", "trip_id"), names(dt_stratum))
  use_cols <- unique(c(id_col, id_metadata, target_raw, predictors_final))
  df <- dt_stratum[, use_cols, with = FALSE] |> as.data.frame()
  
  df[[id_col]]     <- as.factor(df[[id_col]])
  df[[target_raw]] <- suppressWarnings(as.numeric(df[[target_raw]]))
  df$row_id <- seq_len(nrow(df))
  
  if ("bl_hr_person" %in% names(df)) {
    df$bl_hr_person <- suppressWarnings(as.numeric(df$bl_hr_person))
  }
  
  df <- df[is.finite(df[[target_raw]]), , drop = FALSE]
  
  if (nlevels(df[[id_col]]) < V_OUTER) {
    stop(stratum_name, ": too few subjects for outer CV (need >= ", V_OUTER, ").")
  }
  if (!("bl_hr_person" %in% names(df))) {
    stop(stratum_name, ": bl_hr_person missing; baselines require it.")
  }
  
  log_msg(
    stratum_name, ": rows=", nrow(df),
    " | subjects=", nlevels(df[[id_col]]),
    " | predictors=", length(predictors_final)
  )
  log_msg(
    stratum_name, ": Nested grouped CV | outer folds=", V_OUTER,
    " | inner folds=", V_INNER,
    " | participant grouping=", id_col
  )
  
  # Distribution diagnostics
  nhr_tmp <- df[[target_raw]] - df$bl_hr_person
  dist_diag <- tibble(
    stratum = stratum_name,
    n_rows = nrow(df),
    n_subj = nlevels(df[[id_col]]),
    pct_raw_above_baseline = mean(df[[target_raw]] > df$bl_hr_person, na.rm = TRUE),
    mean_nhr = mean(nhr_tmp, na.rm = TRUE),
    sd_nhr   = sd(nhr_tmp, na.rm = TRUE),
    p05_nhr  = as.numeric(quantile(nhr_tmp, 0.05, na.rm = TRUE)),
    p50_nhr  = as.numeric(quantile(nhr_tmp, 0.50, na.rm = TRUE)),
    p95_nhr  = as.numeric(quantile(nhr_tmp, 0.95, na.rm = TRUE))
  )
  
  # Outer CV folds grouped by subject. These folds are used only for final
  # performance estimation. ENet hyperparameters are tuned inside each outer
  # training set using a separate inner grouped CV.
  set.seed(ANALYSIS_SEED)
  outer_folds <- rsample::group_vfold_cv(df, group = !!sym(id_col), v = V_OUTER)
  
  # Fold-safe baselines and context offsets
  offsets_df <- purrr::map2_dfr(
    outer_folds$splits,
    outer_folds$id,
    function(spl, fold_id) {
      tr <- rsample::analysis(spl)
      c_off <- mean(tr[[target_raw]] - tr$bl_hr_person, na.rm = TRUE)
      
      tibble(
        stratum = stratum_name,
        fold = fold_id,
        c_off = c_off,
        n_tr = nrow(tr),
        n_te = nrow(rsample::assessment(spl)),
        n_subj_tr = nlevels(droplevels(tr[[id_col]])),
        n_subj_te = nlevels(droplevels(rsample::assessment(spl)[[id_col]]))
      )
    }
  )
  
  preds_baseline <- purrr::map2_dfr(
    outer_folds$splits,
    outer_folds$id,
    function(spl, fold_id) {
      tr <- rsample::analysis(spl)
      te <- rsample::assessment(spl)
      
      c_off <- mean(tr[[target_raw]] - tr$bl_hr_person, na.rm = TRUE)
      
      bind_rows(
        tibble(
          model = "baseline_only",
          stratum = stratum_name,
          fold = fold_id,
          row_id = te$row_id,
          p_id = as.character(te[[id_col]]),
          dt_time = if ("dt_time" %in% names(te)) te$dt_time else as.POSIXct(NA),
          day_key = if ("day_key" %in% names(te)) te$day_key else NA_integer_,
          trip_id = if ("trip_id" %in% names(te)) as.character(te$trip_id) else NA_character_,
          bl_hr_person = te$bl_hr_person,
          raw_hr_obs = te[[target_raw]],
          raw_hr_hat = te$bl_hr_person
        ),
        tibble(
          model = "baseline_plus_offset",
          stratum = stratum_name,
          fold = fold_id,
          row_id = te$row_id,
          p_id = as.character(te[[id_col]]),
          dt_time = if ("dt_time" %in% names(te)) te$dt_time else as.POSIXct(NA),
          day_key = if ("day_key" %in% names(te)) te$day_key else NA_integer_,
          trip_id = if ("trip_id" %in% names(te)) as.character(te$trip_id) else NA_character_,
          bl_hr_person = te$bl_hr_person,
          raw_hr_obs = te[[target_raw]],
          raw_hr_hat = te$bl_hr_person + c_off
        )
      ) %>%
        mutate(
          nhr_obs = raw_hr_obs - bl_hr_person,
          nhr_hat = raw_hr_hat - bl_hr_person
        )
    }
  )
  
  # Recipe. This recipe is prepped only inside resampling/training data.
  recipe_id_cols <- c(id_col, "row_id", setdiff(id_metadata, if (ALLOW_TRIP_ID_IF_REASONABLE) "trip_id" else character(0)))
  rec <- recipe(raw_hr ~ ., data = df) %>%
    update_role(all_of(recipe_id_cols), new_role = "id") %>%
    step_string2factor(all_nominal_predictors()) %>%
    step_unknown(all_nominal_predictors(), new_level = "Unknown") %>%
    step_novel(all_nominal_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_impute_mode(all_nominal_predictors()) %>%
    step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
    step_zv(all_predictors()) %>%
    step_normalize(all_numeric_predictors())
  
  # Leak check on representative split
  prepped0 <- prep(rec, training = rsample::analysis(outer_folds$splits[[1]]), verbose = FALSE)
  assert_no_forbidden_in_baked(
    prepped0,
    te_df = rsample::assessment(outer_folds$splits[[1]]),
    outcome_name = "raw_hr",
    context = paste0(stratum_name, "_baked_check_outer1")
  )
  
  # ENet spec with tunable penalty and mixture. Tuning happens only in inner CV.
  enet_spec <- linear_reg(penalty = tune(), mixture = tune()) %>% set_engine("glmnet")
  enet_wf   <- workflow() %>% add_recipe(rec) %>% add_model(enet_spec)
  metrics_raw <- metric_set(rmse, mae, rsq)
  
  ctrl_inner <- control_grid(
    save_pred = FALSE,
    save_workflow = FALSE,
    parallel_over = "resamples",
    verbose = FALSE
  )
  
  # Stage 1 grid is fixed in advance. Stage 2 refines within each outer train set.
  grid1 <- tidyr::crossing(
    penalty = 10^seq(-6, -1, length.out = 12),
    mixture = seq(0, 1, length.out = 6)
  )
  
  nested_best_list <- list()
  nested_metrics_list <- list()
  preds_enet_list <- list()
  
  for (k in seq_along(outer_folds$splits)) {
    fold_id <- outer_folds$id[[k]]
    spl <- outer_folds$splits[[k]]
    tr_outer <- rsample::analysis(spl)
    te_outer <- rsample::assessment(spl)
    
    n_subj_outer_train <- nlevels(droplevels(tr_outer[[id_col]]))
    v_inner_use <- min(V_INNER, n_subj_outer_train)
    if (v_inner_use < 2) {
      stop(stratum_name, ": too few outer-training subjects for inner CV in ", fold_id)
    }
    
    set.seed(ANALYSIS_SEED + k)
    inner_folds <- rsample::group_vfold_cv(tr_outer, group = !!sym(id_col), v = v_inner_use)
    
    log_msg(
      stratum_name, ": ", fold_id,
      " | inner tuning Stage1 grid=", nrow(grid1),
      " | inner folds=", v_inner_use,
      " | outer train subjects=", n_subj_outer_train,
      " | outer test subjects=", nlevels(droplevels(te_outer[[id_col]]))
    )
    
    set.seed(ANALYSIS_SEED + k)
    enet_res1 <- tune_grid(
      enet_wf,
      resamples = inner_folds,
      grid = grid1,
      metrics = metrics_raw,
      control = ctrl_inner
    )
    
    best1 <- select_best(enet_res1, metric = "rmse")
    pen <- best1$penalty
    mix <- best1$mixture
    pen_log <- log10(pen)
    pen_grid2 <- 10^seq(pen_log - 0.75, pen_log + 0.75, length.out = 12)
    mix_grid2 <- sort(unique(pmin(1, pmax(0, mix + seq(-0.3, 0.3, length.out = 7)))))
    grid2 <- tidyr::crossing(
      penalty = pen_grid2,
      mixture = mix_grid2
    )
    
    log_msg(
      stratum_name, ": ", fold_id,
      " | Stage1 best penalty=", signif(best1$penalty, 6),
      " mixture=", signif(best1$mixture, 6),
      " | Stage2 grid=", nrow(grid2)
    )
    
    set.seed(ANALYSIS_SEED + 100L + k)
    enet_res2 <- tune_grid(
      enet_wf,
      resamples = inner_folds,
      grid = grid2,
      metrics = metrics_raw,
      control = ctrl_inner
    )
    
    best2 <- select_best(enet_res2, metric = "rmse")
    log_msg(
      stratum_name, ": ", fold_id,
      " | Stage2 best penalty=", signif(best2$penalty, 6),
      " mixture=", signif(best2$mixture, 6)
    )
    
    # Save tuning objects for auditability. These objects do not contain final
    # outer-test performance and are used only to document inner tuning.
    saveRDS(
      list(stage1 = enet_res1, stage2 = enet_res2),
      file.path(out_dir, paste0("enet_nested_tuning_", stratum_name, "_", fold_id, ".rds"))
    )
    
    nested_best_list[[k]] <- tibble(
      stratum = stratum_name,
      fold = fold_id,
      outer_train_rows = nrow(tr_outer),
      outer_test_rows = nrow(te_outer),
      outer_train_subjects = n_subj_outer_train,
      outer_test_subjects = nlevels(droplevels(te_outer[[id_col]])),
      inner_folds = v_inner_use,
      stage1_penalty = best1$penalty,
      stage1_mixture = best1$mixture,
      penalty = best2$penalty,
      mixture = best2$mixture
    )
    
    final_wf_outer <- finalize_workflow(enet_wf, best2)
    final_fit_outer <- fit(final_wf_outer, data = tr_outer)
    
    pred_te <- predict(final_fit_outer, new_data = te_outer) %>%
      dplyr::bind_cols(
        te_outer %>%
          dplyr::transmute(
            row_id,
            p_id = as.character(.data[[id_col]]),
            dt_time = if ("dt_time" %in% names(te_outer)) dt_time else as.POSIXct(NA),
            day_key = if ("day_key" %in% names(te_outer)) day_key else NA_integer_,
            trip_id = if ("trip_id" %in% names(te_outer)) as.character(trip_id) else NA_character_,
            bl_hr_person,
            raw_hr_obs = .data[[target_raw]]
          )
      ) %>%
      dplyr::mutate(
        model = "enet",
        stratum = stratum_name,
        fold = fold_id,
        raw_hr_hat = .pred,
        nhr_obs = raw_hr_obs - bl_hr_person,
        nhr_hat = raw_hr_hat - bl_hr_person
      ) %>%
      dplyr::select(
        model, stratum, fold, row_id, p_id, dt_time, day_key, trip_id,
        bl_hr_person, raw_hr_obs, raw_hr_hat, nhr_obs, nhr_hat
      )
    
    preds_enet_list[[k]] <- pred_te
    
    mm <- compute_metrics_vec(pred_te$raw_hr_obs, pred_te$raw_hr_hat)
    nested_metrics_list[[k]] <- tibble(
      stratum = stratum_name,
      fold = fold_id,
      .metric = c("rmse", "mae", "rsq_cor", "rsq_oos"),
      .estimate = c(mm$rmse, mm$mae, mm$rsq_cor, mm$rsq_oos),
      n = mm$n
    )
  }
  
  best_df <- bind_rows(nested_best_list)
  write_csv(best_df, file.path(out_dir, paste0("best_params_", stratum_name, ".csv")))
  
  enet_cv_metrics <- bind_rows(nested_metrics_list)
  write_csv(enet_cv_metrics, file.path(out_dir, paste0("enet_cv_metrics_", stratum_name, ".csv")))
  
  preds_enet <- bind_rows(preds_enet_list)
  
  # ----------------------------
  # Feature importance: refit full stratum
  # ----------------------------
  # For descriptive coefficient/importance analyses and downstream scenario
  # modeling, refit on the complete stratum using the component-wise median of
  # the hyperparameters selected independently in the outer folds. This full-data
  # refit never contributes to performance estimates, which are based exclusively
  # on held-out outer-fold predictions.
  best_full <- tibble(
    penalty = median(best_df$penalty, na.rm = TRUE),
    mixture = median(best_df$mixture, na.rm = TRUE)
  )
  write_csv(
    tibble(stratum = stratum_name, penalty = best_full$penalty, mixture = best_full$mixture),
    file.path(out_dir, paste0("best_params_full_refit_for_importance_", stratum_name, ".csv"))
  )
  
  log_msg(
    stratum_name, ": refit ENet on FULL data for descriptive importance and downstream analyses | penalty=",
    signif(best_full$penalty, 6), " mixture=", signif(best_full$mixture, 6)
  )
  final_wf_full <- finalize_workflow(enet_wf, best_full)
  final_fit_full <- fit(final_wf_full, data = df)
  saveRDS(
    final_fit_full,
    file.path(out_dir, paste0("enet_full_refit_", stratum_name, ".rds"))
  )
  write_csv(
    tibble(column = names(df), class = map_chr(df, ~paste(class(.x), collapse = "|"))),
    file.path(out_dir, paste0("model_input_schema_", stratum_name, ".csv"))
  )
  
  imp <- extract_enet_importance(
    fitted_wf = final_fit_full,
    penalty_value = best_full$penalty,
    predictors_final = predictors_final
  )
  
  write_csv(imp$term_tbl,  file.path(out_dir, paste0("feature_importance_terms_", stratum_name, ".csv")))
  write_csv(imp$group_tbl, file.path(out_dir, paste0("feature_importance_grouped_", stratum_name, ".csv")))
  
  # Preserve complete coefficient tables for auditability, but create
  # separate presentation tables for non-physiological modulators only.
  # This changes no fitted model or performance estimate.
  imp_terms_modulators <- imp$term_tbl %>%
    filter(group != "bl_hr_person")
  imp_group_modulators <- imp$group_tbl %>%
    filter(group != "bl_hr_person")
  
  write_csv(
    imp_terms_modulators,
    file.path(out_dir, paste0("feature_importance_terms_modulators_only_", stratum_name, ".csv"))
  )
  write_csv(
    imp_group_modulators,
    file.path(out_dir, paste0("feature_importance_grouped_modulators_only_", stratum_name, ".csv"))
  )
  
  p_imp_terms <- plot_importance_bar(
    imp_terms_modulators, name_col = "term", top_n = 30,
    title = paste0("ENet modulator importance (term-level) — ", stratum_name)
  )
  safe_save_pdf(p_imp_terms, file.path(fig_dir, paste0("Fig_FeatureImportance_Terms_", stratum_name, ".pdf")), w = 10.0, h = 7.0)
  safe_save_png(p_imp_terms, file.path(fig_dir, paste0("Fig_FeatureImportance_Terms_", stratum_name, ".png")), w = 10.0, h = 7.0, dpi = 300)
  
  p_imp_group <- plot_importance_bar(
    imp_group_modulators, name_col = "group", top_n = 30,
    title = paste0("ENet modulator importance (grouped) — ", stratum_name)
  )
  safe_save_pdf(p_imp_group, file.path(fig_dir, paste0("Fig_FeatureImportance_Grouped_", stratum_name, ".pdf")), w = 10.0, h = 7.0)
  safe_save_png(p_imp_group, file.path(fig_dir, paste0("Fig_FeatureImportance_Grouped_", stratum_name, ".png")), w = 10.0, h = 7.0, dpi = 300)
  
  # Combine predictions
  preds_df <- bind_rows(preds_baseline, preds_enet)
  
  # Metrics
  m_raw_byfold <- preds_df %>%
    group_by(stratum, model, fold) %>%
    group_modify(~{
      mm <- compute_metrics_vec(.x$raw_hr_obs, .x$raw_hr_hat)
      metrics_to_long(mm)
    }) %>%
    ungroup()
  
  m_raw_overall <- preds_df %>%
    group_by(stratum, model) %>%
    group_modify(~{
      mm <- compute_metrics_vec(.x$raw_hr_obs, .x$raw_hr_hat)
      metrics_to_long(mm)
    }) %>%
    ungroup()
  
  m_nhr_byfold <- preds_df %>%
    group_by(stratum, model, fold) %>%
    group_modify(~{
      mm <- compute_metrics_vec(.x$nhr_obs, .x$nhr_hat)
      metrics_to_long(mm)
    }) %>%
    ungroup()
  
  m_nhr_overall <- preds_df %>%
    group_by(stratum, model) %>%
    group_modify(~{
      mm <- compute_metrics_vec(.x$nhr_obs, .x$nhr_hat)
      metrics_to_long(mm)
    }) %>%
    ungroup()
  
  log_msg("---- Stratum END: ", stratum_name, " ----")
  
  list(
    dist_diag = dist_diag,
    offsets = offsets_df,
    best = best_df,
    preds = preds_df,
    m_raw_byfold = m_raw_byfold,
    m_raw_overall = m_raw_overall,
    m_nhr_byfold = m_nhr_byfold,
    m_nhr_overall = m_nhr_overall
  )
}

# ============================================================
# Run both strata
# ============================================================
res_drive <- run_stratum(dt_drive, "DRIVING", preds_drive)
res_nond  <- run_stratum(dt_nond,  "NONDRIVING_SEDENTARY", preds_nond)

# ============================================================
# Save outputs
# ============================================================
dist_all <- bind_rows(res_drive$dist_diag, res_nond$dist_diag)
write_csv(dist_all, file.path(out_dir, "distribution_diagnostics_by_stratum.csv"))

offsets_all <- bind_rows(res_drive$offsets, res_nond$offsets)
write_csv(offsets_all, file.path(out_dir, "outer_fold_offsets_by_stratum.csv"))

best_all <- bind_rows(res_drive$best, res_nond$best)
write_csv(best_all, file.path(out_dir, "best_params_by_stratum.csv"))

preds_all <- bind_rows(res_drive$preds, res_nond$preds)
write_csv(preds_all, file.path(out_dir, "predictions_all_models_both_strata.csv"))

m_raw_byfold_all  <- bind_rows(res_drive$m_raw_byfold,  res_nond$m_raw_byfold)
m_raw_overall_all <- bind_rows(res_drive$m_raw_overall, res_nond$m_raw_overall)
m_nhr_byfold_all  <- bind_rows(res_drive$m_nhr_byfold,  res_nond$m_nhr_byfold)
m_nhr_overall_all <- bind_rows(res_drive$m_nhr_overall, res_nond$m_nhr_overall)

write_csv(
  m_raw_byfold_all,
  file.path(out_dir, "compare_metrics_rawhr_by_fold_by_stratum.csv")
)

# Pooled out-of-fold metrics weight observations rather than folds. They are
# useful diagnostics but are not the summaries shown in Figure 5 or Table 4.
write_csv(
  m_raw_overall_all,
  file.path(out_dir, "compare_metrics_rawhr_pooled_oof_by_stratum.csv")
)
write_csv(
  m_nhr_byfold_all,
  file.path(out_dir, "compare_metrics_nhrhat_by_fold_by_stratum.csv")
)
write_csv(
  m_nhr_overall_all,
  file.path(out_dir, "compare_metrics_nhrhat_pooled_oof_by_stratum.csv")
)

# Manuscript Figure 5 and Table 4 use the unweighted mean and standard
# error of the five held-out outer-fold metrics.
manuscript_metrics <- m_raw_byfold_all %>%
  group_by(stratum, model, .metric) %>%
  summarise(
    mean = mean(.estimate, na.rm = TRUE),
    sd = sd(.estimate, na.rm = TRUE),
    se = sd / sqrt(sum(is.finite(.estimate))),
    n_outer_folds = sum(is.finite(.estimate)),
    .groups = "drop"
  )

write_csv(
  manuscript_metrics,
  file.path(out_dir, "manuscript_table4_metrics.csv")
)

log_msg("Saved result CSVs, including the fold summaries used in manuscript Table 4.")

# ============================================================
# Figures (performance)
# ============================================================
p_above <- dist_all %>%
  ggplot(aes(x = stratum, y = pct_raw_above_baseline)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "% RAW_HR > baseline",
    title = "How often is HR above baseline? (Driving vs Non-driving sedentary)"
  ) +
  theme_minimal(base_size = 12)
safe_save_pdf(p_above, file.path(fig_dir, "Fig_AboveBaselineRate_byStratum.pdf"), w = 7.4, h = 4.2)
safe_save_png(p_above, file.path(fig_dir, "Fig_AboveBaselineRate_byStratum.png"), w = 7.4, h = 4.2, dpi = 300)

p_perf_raw <- manuscript_metrics %>%
  ggplot(aes(x = model, y = mean)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  facet_grid(.metric ~ stratum, scales = "free_y") +
  labs(
    x = NULL,
    y = "Mean across outer folds",
    title = "Nested-CV performance — raw heart rate",
    subtitle = "Error bars show standard errors across five outer folds"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
safe_save_pdf(
  p_perf_raw,
  file.path(fig_dir, "Fig_Performance_ByOuterFold_RAWHR_byStratum.pdf"),
  w = 10.0, h = 5.2
)
safe_save_png(
  p_perf_raw,
  file.path(fig_dir, "Fig_Performance_ByOuterFold_RAWHR_byStratum.png"),
  w = 10.0, h = 5.2, dpi = 300
)

manuscript_metrics_nhr <- m_nhr_byfold_all %>%
  group_by(stratum, model, .metric) %>%
  summarise(
    mean = mean(.estimate, na.rm = TRUE),
    sd = sd(.estimate, na.rm = TRUE),
    se = sd / sqrt(sum(is.finite(.estimate))),
    n_outer_folds = sum(is.finite(.estimate)),
    .groups = "drop"
  )

p_perf_nhr <- manuscript_metrics_nhr %>%
  ggplot(aes(x = model, y = mean)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  facet_grid(.metric ~ stratum, scales = "free_y") +
  labs(
    x = NULL,
    y = "Mean across outer folds",
    title = "Nested-CV performance — baseline-referenced heart rate",
    subtitle = "Error bars show standard errors across five outer folds"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
safe_save_pdf(
  p_perf_nhr,
  file.path(fig_dir, "Fig_Performance_ByOuterFold_NHR_byStratum.pdf"),
  w = 10.0, h = 5.2
)
safe_save_png(
  p_perf_nhr,
  file.path(fig_dir, "Fig_Performance_ByOuterFold_NHR_byStratum.png"),
  w = 10.0, h = 5.2, dpi = 300
)

# Compact output manifest for repository users and downstream scripts.
# It is written last and therefore lists all completed outputs except itself.
manifest <- tibble(
  file = list.files(out_dir, recursive = TRUE),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
write_csv(manifest, file.path(out_dir, "output_manifest.csv"))

log_msg("DONE. Outputs in: ", out_dir)