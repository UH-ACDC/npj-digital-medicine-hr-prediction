# ============================================================
# 07_figure7_table5_long_horizon_scenarios.R
#
# PURPOSE
#   Reproduce Figure 7 and Table 5 for:
#
#     Wearable sensing reveals the structure of cardiac
#     activation associated with everyday driving
#
#   Figure 7A:
#     Illustrative trait-anxiety scenario projections of annual
#     cumulative baseline-referenced cardiac activation.
#
#   Figure 7B:
#     Normalized out-of-fold ENet RMSE after temporal aggregation
#     over contiguous 1-, 2-, 5-, 10-, 15-, 30-, and 60-min windows.
#
#   Table 5:
#     Minute-level predicted activation and annual cumulative
#     NHR-hours for the lowest and highest observed trait-anxiety
#     values, with participant-bootstrap 5th-95th percentile ranges.
#
# MANUSCRIPT ALIGNMENT
#   - Primary resolution is fixed at 60 s.
#   - The final DRIVING ENet is refit on the complete DRIVING
#     stratum using the manuscript-selected hyperparameters from
#     best_params_DRIVING.csv.
#   - Every observed DRIVING epoch is evaluated twice, with trait
#     anxiety fixed to the lowest and highest observed participant
#     values; all other predictors and participant baselines remain
#     unchanged.
#   - Uncertainty is obtained by paired participant bootstrap.
#   - Figure 7B uses held-out outer-fold ENet predictions only.
#
# INPUTS
#   Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
#   A compatible run folder under Results/nubi_ml containing:
#     predictor_list_DRIVING.csv
#     best_params_DRIVING.csv
#     best_params_DRIVING.csv
#     one compatible all-model out-of-fold prediction CSV
#
# OUTPUTS
#   Results/paper_figs/<timestamp>_60sec_figure7_table5_long_horizon/
#
# REPOSITORY SCOPE
#   The public repository begins from the final clean MASTER dataset
#   and curated nested-CV outputs. It does not reconstruct raw sensor
#   ingestion or the full upstream model-selection pipeline.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(tidymodels)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

options(warn = 1)
if (exists("tidymodels_prefer", mode = "function")) tidymodels_prefer()
set.seed(20260309)

# ============================================================
# USER SETTINGS
# ============================================================

LOCAL_TZ <- "America/Chicago"
RES_SECONDS <- 60L

# Use an exact folder name inside Results/nubi_ml, or leave NULL
# to select the most recently modified compatible 60-s run.
RUN_DIR_NAME <- NULL

TRAIT_VAR_CANDIDATES <- c(
  "trait_anxiety",
  "stai_trait",
  "traitanxiety",
  "stai_trait_total"
)

TRAIT_LABEL_LOW  <- "Lowest trait anxiety"
TRAIT_LABEL_HIGH <- "Highest trait anxiety"

DAYS_PER_YEAR <- 365
SCENARIOS_HOURS_PER_DAY <- c(
  "30 min/day commute"            = 0.5,
  "2 hr/day commute"              = 2.0,
  "8 hr/day professional driving" = 8.0
)

N_BOOT <- 20000L
HORIZONS_MIN <- c(1L, 2L, 5L, 10L, 15L, 30L, 60L)
MAX_GAP_MULTIPLIER <- 1.5

PAL_TRAIT <- c(
  "Lowest trait anxiety"  = "#4DD3D3",
  "Highest trait anxiety" = "#F4A3A3"
)

PAL_STRATUM <- c(
  "Driving" = "#E69F00",
  "Non-driving sedentary" = "gray45"
)

PNG_DPI <- 300

# Manuscript reference values used only for an end-of-run
# reproducibility check. The script does not force these values.
MANUSCRIPT_LOW_MINUTE_NHR  <- 10.91
MANUSCRIPT_HIGH_MINUTE_NHR <- 11.94
MANUSCRIPT_DIFF_MINUTE_NHR <- 1.03
MANUSCRIPT_CHECK_TOL_BPM   <- 0.02

# ============================================================
# PATHS
# ============================================================

this_script <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) NA_character_
)

script_dir <- if (!is.na(this_script) && file.exists(this_script)) {
  dirname(this_script)
} else {
  getwd()
}
script_dir <- normalizePath(script_dir, mustWork = TRUE)
setwd(script_dir)

project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_path <- file.path(
  project_root,
  "Data",
  sprintf("NUBI_Data_%dsec_Level_MASTER_CLEAN.csv", RES_SECONDS)
)
ml_root <- file.path(project_root, "Results", "nubi_ml")
paper_fig_root <- file.path(project_root, "Results", "paper_figs")

if (!file.exists(data_path)) stop("Missing MASTER dataset: ", data_path)
if (!dir.exists(ml_root)) stop("Missing ML results directory: ", ml_root)
dir.create(paper_fig_root, recursive = TRUE, showWarnings = FALSE)

message("Script directory: ", script_dir)
message("Project root: ", project_root)
message("Resolution: ", RES_SECONDS, " sec")

# ============================================================
# HELPERS
# ============================================================

snakeify <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

norm_chr <- function(x) trimws(tolower(as.character(x)))

parse_time_local <- function(x, tz = LOCAL_TZ) {
  if (inherits(x, "POSIXt")) return(with_tz(x, tzone = tz))
  x <- as.character(x)
  z <- suppressWarnings(ymd_hms(x, tz = tz, quiet = TRUE))
  if (!all(is.na(z))) return(z)
  suppressWarnings(parse_date_time(
    x,
    orders = c("ymd HMS", "ymd HM", "mdy HMS", "mdy HM",
               "dmy HMS", "dmy HM"),
    tz = tz
  ))
}

canonicalize_names <- function(dt) {
  stopifnot(is.data.table(dt))

  old <- names(dt)
  sn <- snakeify(old)

  aliases <- c(
    pid = "p_id",
    participant_id = "p_id",
    participantid = "p_id",
    timestamp = "time",
    datetime = "time",
    date_time = "time",
    datasource = "activity",
    data_source = "activity",
    source = "activity",
    hr_bl = "bl_hr",
    hrbl = "bl_hr",
    baseline = "bl_hr",
    daynum = "day_num",
    weather = "weather_info",
    traitanxiety = "trait_anxiety"
  )

  new <- sn
  hit <- sn %in% names(aliases)
  new[hit] <- unname(aliases[sn[hit]])

  if (!identical(old, new)) {
    new <- make.unique(new, sep = "_")
    setnames(dt, old, new)
  }

  dt
}

normalize_stratum <- function(x) {
  z <- toupper(trimws(as.character(x)))
  z <- str_replace_all(z, "[[:space:]/-]+", "_")

  case_when(
    z %in% c("DRIVING", "DRIVE") ~ "DRIVING",
    z %in% c(
      "NONDRIVING_SEDENTARY",
      "NON_DRIVING_SEDENTARY",
      "NONDRIVING",
      "NON_DRIVING",
      "SEDENTARY_NONDRIVING",
      "SEDENTARY_NON_DRIVING",
      "NONDRIVINGSEDENTARY"
    ) ~ "NONDRIVING_SEDENTARY",
    TRUE ~ NA_character_
  )
}

normalize_model <- function(x) {
  z <- tolower(trimws(as.character(x)))
  z <- str_replace_all(z, "[[:space:]/-]+", "_")

  case_when(
    z %in% c("enet", "elastic_net", "elasticnet", "full_enet") ~ "enet",
    z %in% c(
      "baseline_offset", "baseline_plus_offset",
      "baseline_context_offset", "offset"
    ) ~ "baseline_offset",
    z %in% c(
      "baseline0", "baseline_0", "baseline", "baseline_only",
      "baselineonly"
    ) ~ "baseline0",
    TRUE ~ NA_character_
  )
}

first_existing_col <- function(df, candidates, required = TRUE, what = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0L) {
    if (required) {
      stop("Could not find ", what, ". Tried: ",
           paste(candidates, collapse = ", "))
    }
    return(NA_character_)
  }
  hit[1]
}

safe_save_pdf <- function(plot_obj, path, w, h) {
  ok <- tryCatch({
    ggsave(
      filename = path,
      plot = plot_obj,
      width = w,
      height = h,
      device = grDevices::cairo_pdf
    )
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    ggsave(
      filename = path,
      plot = plot_obj,
      width = w,
      height = h,
      device = "pdf",
      useDingbats = FALSE
    )
  }
}

find_exact_file <- function(root, basename_required) {
  hits <- list.files(
    root,
    pattern = paste0("^", basename_required, "$"),
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(hits) == 0L) return(NA_character_)

  rel <- substring(hits, nchar(root) + 2L)
  depth <- lengths(strsplit(rel, .Platform$file.sep, fixed = TRUE))
  info <- file.info(hits)
  hits[order(depth, -as.numeric(info$mtime))[1]]
}

prediction_basenames <- c(
  "predictions_all_models_both_strata.csv",
  "predictions_all_models.csv",
  "oof_predictions_all_models_both_strata.csv",
  "oof_predictions_all_models.csv"
)

resolve_run_inputs <- function(run_dir) {
  predictor <- find_exact_file(run_dir, "predictor_list_DRIVING.csv")

  # Figure 7 and Table 5 must use the manuscript-selected
  # hyperparameters. The importance-refit file belongs to the
  # Figure 6 interpretation workflow and is not interchangeable.
  param_file <- find_exact_file(run_dir, "best_params_DRIVING.csv")

  # Retain the Figure 6 full-refit parameter file only as a
  # diagnostic comparison, never as a fallback for Figure 7.
  importance_param_file <- find_exact_file(
    run_dir,
    "best_params_full_refit_for_importance_DRIVING.csv"
  )

  pred_hits <- vapply(
    prediction_basenames,
    function(x) find_exact_file(run_dir, x),
    character(1)
  )
  pred_file <- unname(pred_hits[!is.na(pred_hits)][1])
  if (length(pred_file) == 0L) pred_file <- NA_character_

  c(
    predictor = predictor,
    params = param_file,
    importance_params = importance_param_file,
    predictions = pred_file
  )
}

auto_pick_run_dir <- function(root, res_seconds) {
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.info(dirs)$isdir %in% TRUE]
  res_pat <- paste0("(^|[^0-9])", res_seconds, "sec([^0-9]|$)")
  dirs <- dirs[grepl(res_pat, basename(dirs), ignore.case = TRUE, perl = TRUE)]

  if (length(dirs) == 0L) {
    stop("No ", res_seconds, "-s run folders found under: ", root)
  }

  ok <- vapply(
    dirs,
    function(d) {
      p <- resolve_run_inputs(d)
      all(!is.na(p[c("predictor", "params", "predictions")]))
    },
    logical(1)
  )

  cand <- dirs[ok]
  if (length(cand) == 0L) {
    inventory <- vapply(
      dirs,
      function(d) {
        p <- resolve_run_inputs(d)
        paste0(
          basename(d), ": predictor=", !is.na(p["predictor"]),
          ", manuscript_params=", !is.na(p["params"]),
          ", importance_params=", !is.na(p["importance_params"]),
          ", predictions=", !is.na(p["predictions"])
        )
      },
      character(1)
    )
    stop(
      "Could not find a compatible 60-s ML run.\n",
      paste(inventory, collapse = "\n")
    )
  }

  cand[which.max(file.info(cand)$mtime)]
}

make_day_key <- function(dt) {
  if ("day_num" %in% names(dt)) {
    x <- suppressWarnings(as.integer(str_extract(as.character(dt$day_num), "\\d+")))
    if (!all(is.na(x))) return(x)
  }
  if ("days" %in% names(dt)) {
    x <- suppressWarnings(as.integer(str_extract(as.character(dt$days), "\\d+")))
    if (!all(is.na(x))) return(x)
  }
  as.integer(as.Date(dt$dt_time, tz = LOCAL_TZ))
}

add_bl_hr_person <- function(dt) {
  by_day <- dt[
    is.finite(bl_hr) & !is.na(day_key),
    .(bl_hr_day = median(bl_hr, na.rm = TRUE)),
    by = .(p_id, day_key)
  ]

  if (nrow(by_day) == 0L) {
    stop("Could not calculate participant-day baselines.")
  }

  by_person <- by_day[
    ,
    .(bl_hr_person = mean(bl_hr_day, na.rm = TRUE)),
    by = p_id
  ]

  merge(dt, by_person, by = "p_id", all.x = TRUE)
}

dyn_base_vars <- c(
  "speed", "ff", "ff_speed", "atp", "rtp", "jf",
  "energy_acc", "energy_rot"
)

add_dynamics <- function(dt,
                         id_col = "p_id",
                         time_col = "dt_time",
                         vars = dyn_base_vars,
                         res_seconds = RES_SECONDS,
                         windows_min = c(1, 3, 5)) {
  stopifnot(is.data.table(dt))
  vars <- intersect(vars, names(dt))
  if (length(vars) == 0L) return(dt)

  setorderv(dt, c(id_col, time_col))

  for (v in vars) {
    dt[, (v) := suppressWarnings(as.numeric(get(v)))]

    dt[, paste0(v, "_lag1") := shift(get(v), 1L), by = id_col]
    dt[, paste0(v, "_diff1") := get(v) - get(paste0(v, "_lag1")),
       by = id_col]

    for (wm in windows_min) {
      k <- max(2L, as.integer(round((wm * 60) / res_seconds)))
      tag <- paste0("_", wm, "m")

      rm_name <- paste0(v, "_rm", tag)
      rs_name <- paste0(v, "_rs", tag)
      sl_name <- paste0(v, "_slope", tag)

      dt[, (rm_name) := frollmean(
        get(v), n = k, align = "right", fill = NA_real_
      ), by = id_col]

      dt[, (rs_name) := {
        m1 <- frollmean(get(v), n = k, align = "right", fill = NA_real_)
        m2 <- frollmean(get(v)^2, n = k, align = "right", fill = NA_real_)
        sqrt(pmax(m2 - m1^2, 0))
      }, by = id_col]

      dt[, (sl_name) := (
        get(v) - shift(get(v), k - 1L)
      ) / ((k - 1L) * res_seconds), by = id_col]
    }
  }

  dt
}

# ============================================================
# SELECT RUN AND CREATE OUTPUT DIRECTORY
# ============================================================

run_dir <- if (!is.null(RUN_DIR_NAME)) {
  file.path(ml_root, RUN_DIR_NAME)
} else {
  auto_pick_run_dir(ml_root, RES_SECONDS)
}

if (!dir.exists(run_dir)) stop("Selected run folder does not exist: ", run_dir)
run_inputs <- resolve_run_inputs(run_dir)

message("Selected ML run: ", basename(run_dir))
message("Predictor list: ", run_inputs["predictor"])
message("Manuscript parameters: ", run_inputs["params"])
message("Importance-refit parameters (diagnostic only): ", run_inputs["importance_params"])
message("OOF predictions: ", run_inputs["predictions"])

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(
  paper_fig_root,
  paste0(stamp, "_", RES_SECONDS, "sec_figure7_table5_long_horizon")
)
fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste(..., collapse = "")
  )
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Script: 07_figure7_table5_long_horizon_scenarios.R")
log_msg("Project root: ", project_root)
log_msg("Resolution: ", RES_SECONDS, " sec")
log_msg("MASTER data: ", data_path)
log_msg("ML run: ", run_dir)
log_msg("Predictor list: ", run_inputs["predictor"])
log_msg("Manuscript parameter file: ", run_inputs["params"])
log_msg("Importance-refit parameter file (diagnostic only): ", run_inputs["importance_params"])
log_msg("OOF prediction file: ", run_inputs["predictions"])
log_msg("N_BOOT: ", N_BOOT)
log_msg("Output directory: ", out_dir)

# ============================================================
# READ SAVED MODEL SPECIFICATION
# ============================================================

pred_list <- read_csv(run_inputs["predictor"], show_col_types = FALSE)
pred_col <- first_existing_col(
  pred_list,
  c("predictor", "term", "variable", "feature"),
  TRUE,
  "predictor-list column"
)
preds_drive <- unique(as.character(pred_list[[pred_col]]))
preds_drive <- preds_drive[!is.na(preds_drive) & nzchar(preds_drive)]

param_df <- read_csv(run_inputs["params"], show_col_types = FALSE)
if (!all(c("penalty", "mixture") %in% names(param_df))) {
  stop("Parameter file must contain penalty and mixture columns: ",
       run_inputs["params"])
}

# Exact manuscript behavior: the original Figure 7A script used
# the first saved row in best_params_DRIVING.csv.
penalty_val <- suppressWarnings(as.numeric(param_df$penalty[1]))
mixture_val <- suppressWarnings(as.numeric(param_df$mixture[1]))

if (!is.finite(penalty_val) || !is.finite(mixture_val)) {
  stop("Invalid first-row penalty/mixture values in: ", run_inputs["params"])
}

log_msg(
  "Manuscript Figure 7A hyperparameters (first saved row): penalty=",
  penalty_val, " | mixture=", mixture_val
)

write_csv(
  tibble(
    resolution_seconds = RES_SECONDS,
    master_file = basename(data_path),
    run_folder = basename(run_dir),
    predictor_file = basename(run_inputs["predictor"]),
    parameter_file = basename(run_inputs["params"]),
    parameter_row_used = 1L,
    oof_prediction_file = basename(run_inputs["predictions"]),
    penalty = penalty_val,
    mixture = mixture_val,
    bootstrap_replicates = N_BOOT
  ),
  file.path(out_dir, "figure7_model_specification.csv")
)

# ============================================================
# READ MASTER DATA AND RECREATE MODEL FRAME
# ============================================================

# This block intentionally follows the manuscript-generating Figure 7A
# script. In particular, it uses the activity column for the DRIVING
# stratum and includes only saved predictors already present in the clean
# MASTER dataset; it does not reconstruct additional dynamic predictors.
dt <- fread(data_path, showProgress = TRUE)
dt <- canonicalize_names(dt)

required_master <- c("p_id", "time", "raw_hr", "activity", "bl_hr")
missing_master <- setdiff(required_master, names(dt))
if (length(missing_master) > 0L) {
  stop("MASTER data missing: ", paste(missing_master, collapse = ", "))
}

if (!("dt_time" %in% names(dt)) || !inherits(dt$dt_time, "POSIXt")) {
  dt[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]
}
dt <- dt[!is.na(dt_time)]

dt[, activity_norm := norm_chr(activity)]
dt[, raw_hr_num := suppressWarnings(as.numeric(raw_hr))]
dt[, bl_hr_num := suppressWarnings(as.numeric(bl_hr))]
dt <- dt[is.finite(raw_hr_num)]

dt[, day_key := make_day_key(dt)]

baseline_by_subj_day <- dt[
  is.finite(bl_hr_num) & !is.na(day_key),
  .(bl_hr_day = median(bl_hr_num, na.rm = TRUE)),
  by = .(p_id, day_key)
]
if (nrow(baseline_by_subj_day) == 0L) {
  stop("Could not calculate participant-day baselines.")
}

baseline_by_subj <- baseline_by_subj_day[
  , .(bl_hr_person = mean(bl_hr_day, na.rm = TRUE)),
  by = p_id
]

dt <- merge(dt, baseline_by_subj, by = "p_id", all.x = TRUE)
dt <- dt[is.finite(bl_hr_person)]

dt_drive <- copy(dt[activity_norm == "driving"])
if (nrow(dt_drive) == 0L) stop("No DRIVING rows found.")

preds_drive_use <- intersect(preds_drive, names(dt_drive))
missing_saved_predictors <- setdiff(preds_drive, preds_drive_use)
if (length(missing_saved_predictors) > 0L) {
  log_msg(
    "Saved predictors absent from MASTER and omitted exactly as in the original Figure 7A script: ",
    paste(missing_saved_predictors, collapse = ", ")
  )
}

trait_candidates_present <- intersect(TRAIT_VAR_CANDIDATES, names(dt_drive))
model_cols <- unique(c(
  "p_id", "raw_hr_num", "bl_hr_person", "weather_info",
  trait_candidates_present, preds_drive_use, "day_key"
))
model_cols <- intersect(model_cols, names(dt_drive))

model_df <- as.data.frame(dt_drive[, ..model_cols])
model_df$p_id <- factor(model_df$p_id)
model_df$raw_hr <- model_df$raw_hr_num
model_df$raw_hr_num <- NULL

trait_var <- TRAIT_VAR_CANDIDATES[TRAIT_VAR_CANDIDATES %in% names(model_df)][1]
if (length(trait_var) == 0L || is.na(trait_var)) {
  stop(
    "Could not find a supported trait-anxiety variable. Tried: ",
    paste(TRAIT_VAR_CANDIDATES, collapse = ", ")
  )
}

log_msg(
  "Full-refit DRIVING rows: ", nrow(model_df),
  " | participants: ", n_distinct(model_df$p_id),
  " | saved predictors used: ", length(preds_drive_use)
)
log_msg("Trait variable: ", trait_var)

trait_vals <- suppressWarnings(as.numeric(model_df[[trait_var]]))
if (all(!is.finite(trait_vals))) {
  stop("Trait variable has no finite numeric values: ", trait_var)
}
model_df[[trait_var]] <- trait_vals

trait_audit <- model_df %>%
  filter(is.finite(.data[[trait_var]])) %>%
  group_by(p_id) %>%
  summarise(
    trait_min = min(.data[[trait_var]], na.rm = TRUE),
    trait_max = max(.data[[trait_var]], na.rm = TRUE),
    trait_range = trait_max - trait_min,
    n_rows = n(),
    .groups = "drop"
  )

write_csv(
  trait_audit,
  file.path(out_dir, "figure7_trait_value_consistency_by_subject.csv")
)

if (any(!is.finite(trait_audit$trait_range)) ||
    any(trait_audit$trait_range > 1e-8, na.rm = TRUE)) {
  stop("Trait anxiety is not constant within participant.")
}

trait_subject <- model_df %>%
  filter(is.finite(.data[[trait_var]])) %>%
  arrange(p_id) %>%
  group_by(p_id) %>%
  summarise(
    trait_value = dplyr::first(.data[[trait_var]]),
    bl_hr_person = dplyr::first(bl_hr_person),
    n_rows_drive = n(),
    n_days_drive = n_distinct(day_key),
    .groups = "drop"
  ) %>%
  arrange(trait_value)

trait_low <- trait_subject$trait_value[1]
trait_high <- trait_subject$trait_value[nrow(trait_subject)]

log_msg("Observed trait-anxiety extremes: ", trait_low, " and ", trait_high)

# ============================================================
# REFIT FINAL DRIVING ENET
# ============================================================

rec <- recipe(raw_hr ~ ., data = model_df) %>%
  update_role(p_id, new_role = "id") %>%
  step_string2factor(all_nominal_predictors()) %>%
  step_unknown(all_nominal_predictors(), new_level = "Unknown") %>%
  step_novel(all_nominal_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors(), one_hot = TRUE) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_numeric_predictors())

enet_spec <- linear_reg(
  penalty = penalty_val,
  mixture = mixture_val
) %>%
  set_engine("glmnet")

wf <- workflow() %>%
  add_recipe(rec) %>%
  add_model(enet_spec)

log_msg("Fitting final ENet on complete DRIVING stratum...")
fit_full <- fit(wf, data = model_df)

predict_nhr <- function(fitted_workflow, new_data) {
  nd <- new_data
  pr <- predict(fitted_workflow, new_data = nd)$.pred
  as.numeric(pr - nd$bl_hr_person)
}

make_counterfactual <- function(base_df, trait_value, condition_label) {
  nd <- base_df
  nd[[trait_var]] <- trait_value

  tibble(
    condition = condition_label,
    p_id = as.character(nd$p_id),
    row_id = seq_len(nrow(nd)),
    trait_value = trait_value,
    bl_hr_person = nd$bl_hr_person,
    nhr_hat = predict_nhr(fit_full, nd)
  )
}

cf_rows <- bind_rows(
  make_counterfactual(model_df, trait_low, TRAIT_LABEL_LOW),
  make_counterfactual(model_df, trait_high, TRAIT_LABEL_HIGH)
)

write_csv(
  cf_rows,
  file.path(out_dir, "figure7_trait_counterfactual_predictions_by_row.csv")
)

cf_subject <- cf_rows %>%
  group_by(condition, p_id) %>%
  summarise(
    n_rows = n(),
    mean_nhr_hat = mean(nhr_hat, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  cf_subject,
  file.path(out_dir, "figure7_trait_counterfactual_subject_summary.csv")
)

# Exact manuscript behavior: bootstrap the low- and high-trait
# counterfactual subject summaries separately. The point estimates are
# unaffected, while the interval draws reproduce the original script.
bootstrap_counterfactual_subjects <- function(subject_summary, B) {
  ids <- unique(subject_summary$p_id)
  if (length(ids) < 3L) stop("Too few participants for bootstrap.")

  bind_rows(lapply(seq_len(B), function(b) {
    sampled_ids <- sample(ids, size = length(ids), replace = TRUE)
    dd <- tibble(p_id = sampled_ids) %>%
      left_join(subject_summary, by = "p_id")

    tibble(
      sim = b,
      mean_nhr_hat = mean(dd$mean_nhr_hat, na.rm = TRUE)
    )
  }))
}

boot_long <- bind_rows(
  bootstrap_counterfactual_subjects(
    cf_subject %>% filter(condition == TRAIT_LABEL_LOW),
    N_BOOT
  ) %>% mutate(condition = TRAIT_LABEL_LOW),
  bootstrap_counterfactual_subjects(
    cf_subject %>% filter(condition == TRAIT_LABEL_HIGH),
    N_BOOT
  ) %>% mutate(condition = TRAIT_LABEL_HIGH)
)

write_csv(
  boot_long,
  file.path(out_dir, "figure7_trait_counterfactual_bootstrap_draws.csv")
)

# Separate bootstrap streams were used in the manuscript script. Pair them
# by replicate number only to form the reported high-minus-low summaries.
boot_wide <- boot_long %>%
  select(sim, condition, mean_nhr_hat) %>%
  pivot_wider(names_from = condition, values_from = mean_nhr_hat) %>%
  transmute(
    sim,
    low = .data[[TRAIT_LABEL_LOW]],
    high = .data[[TRAIT_LABEL_HIGH]]
  )

scenarios <- tibble(
  scenario = names(SCENARIOS_HOURS_PER_DAY),
  hours_per_day = as.numeric(SCENARIOS_HOURS_PER_DAY),
  annual_hours = hours_per_day * DAYS_PER_YEAR
)

annual_draws <- crossing(
  boot_long,
  scenarios
) %>%
  mutate(annual_nhr_hours = mean_nhr_hat * annual_hours)

scenario_summary <- annual_draws %>%
  group_by(scenario, hours_per_day, annual_hours, condition) %>%
  summarise(
    minute_mean = mean(mean_nhr_hat),
    minute_q05 = quantile(mean_nhr_hat, 0.05),
    minute_q95 = quantile(mean_nhr_hat, 0.95),
    annual_mean = mean(annual_nhr_hours),
    annual_q05 = quantile(annual_nhr_hours, 0.05),
    annual_q95 = quantile(annual_nhr_hours, 0.95),
    .groups = "drop"
  ) %>%
  mutate(
    scenario = factor(
      scenario,
      levels = names(SCENARIOS_HOURS_PER_DAY)
    ),
    condition = factor(
      condition,
      levels = c(TRAIT_LABEL_LOW, TRAIT_LABEL_HIGH)
    )
  ) %>%
  arrange(scenario, condition)

write_csv(
  scenario_summary,
  file.path(out_dir, "figure7_trait_counterfactual_load_summary.csv")
)

# Paired high-low differences.
difference_draws <- boot_wide %>%
  mutate(minute_difference = high - low)

difference_summary <- scenarios %>%
  mutate(
    scenario = factor(
      scenario,
      levels = names(SCENARIOS_HOURS_PER_DAY)
    ),
    minute_difference = mean(difference_draws$minute_difference),
    annual_difference = minute_difference * annual_hours
  ) %>%
  arrange(scenario)

write_csv(
  difference_summary,
  file.path(out_dir, "figure7_trait_counterfactual_differences.csv")
)

# Table 5 support file, including the minute-level row.
minute_summary <- scenario_summary %>%
  filter(scenario == names(SCENARIOS_HOURS_PER_DAY)[1]) %>%
  select(condition, minute_mean, minute_q05, minute_q95)

table5_minute <- minute_summary %>%
  select(condition, mean = minute_mean, q05 = minute_q05, q95 = minute_q95) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(mean, q05, q95)
  ) %>%
  transmute(
    driving_schedule = "Minute-level predicted activation (bpm)",
    lowest_trait_anxiety = sprintf(
      "%.2f [%.2f-%.2f]",
      .data[[paste0("mean_", TRAIT_LABEL_LOW)]],
      .data[[paste0("q05_", TRAIT_LABEL_LOW)]],
      .data[[paste0("q95_", TRAIT_LABEL_LOW)]]
    ),
    highest_trait_anxiety = sprintf(
      "%.2f [%.2f-%.2f]",
      .data[[paste0("mean_", TRAIT_LABEL_HIGH)]],
      .data[[paste0("q05_", TRAIT_LABEL_HIGH)]],
      .data[[paste0("q95_", TRAIT_LABEL_HIGH)]]
    ),
    high_minus_low = sprintf(
      "%.2f",
      mean(difference_draws$minute_difference)
    )
  )

table5_annual <- scenario_summary %>%
  mutate(
    scenario = factor(
      scenario,
      levels = names(SCENARIOS_HOURS_PER_DAY)
    )
  ) %>%
  arrange(scenario, condition) %>%
  select(
    scenario, condition,
    annual_mean, annual_q05, annual_q95
  ) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(annual_mean, annual_q05, annual_q95)
  ) %>%
  left_join(
    difference_summary %>% select(scenario, annual_difference),
    by = "scenario"
  ) %>%
  transmute(
    driving_schedule = paste0(scenario, " (NHR-hours/year)"),
    lowest_trait_anxiety = sprintf(
      "%.0f [%.0f-%.0f]",
      .data[[paste0("annual_mean_", TRAIT_LABEL_LOW)]],
      .data[[paste0("annual_q05_", TRAIT_LABEL_LOW)]],
      .data[[paste0("annual_q95_", TRAIT_LABEL_LOW)]]
    ),
    highest_trait_anxiety = sprintf(
      "%.0f [%.0f-%.0f]",
      .data[[paste0("annual_mean_", TRAIT_LABEL_HIGH)]],
      .data[[paste0("annual_q05_", TRAIT_LABEL_HIGH)]],
      .data[[paste0("annual_q95_", TRAIT_LABEL_HIGH)]]
    ),
    high_minus_low = sprintf("%.0f", annual_difference)
  )

table5 <- bind_rows(table5_minute, table5_annual)

write_csv(
  table5,
  file.path(out_dir, "table5_long_horizon_scenario_projections.csv")
)

write_csv(
  tibble(
    condition = c(TRAIT_LABEL_LOW, TRAIT_LABEL_HIGH),
    trait_variable = trait_var,
    trait_value = c(trait_low, trait_high),
    n_driving_rows_evaluated = nrow(model_df),
    n_participants = n_distinct(model_df$p_id),
    baseline_handling = "Observed participant-specific baseline retained",
    other_predictors = "All other observed row-level predictors retained"
  ),
  file.path(out_dir, "figure7_trait_extreme_parameters.csv")
)


# ============================================================
# MANUSCRIPT REPRODUCIBILITY CHECK
# ============================================================

minute_check <- scenario_summary %>%
  filter(
    as.character(scenario) == names(SCENARIOS_HOURS_PER_DAY)[1]
  ) %>%
  select(condition, minute_mean) %>%
  mutate(condition = as.character(condition)) %>%
  pivot_wider(
    names_from = condition,
    values_from = minute_mean
  )

regen_low <- minute_check[[TRAIT_LABEL_LOW]]
regen_high <- minute_check[[TRAIT_LABEL_HIGH]]
regen_diff <- mean(difference_draws$minute_difference)

check_tbl <- tibble(
  quantity = c(
    "Lowest trait-anxiety minute NHR",
    "Highest trait-anxiety minute NHR",
    "High-minus-low minute NHR"
  ),
  manuscript_value = c(
    MANUSCRIPT_LOW_MINUTE_NHR,
    MANUSCRIPT_HIGH_MINUTE_NHR,
    MANUSCRIPT_DIFF_MINUTE_NHR
  ),
  regenerated_value = c(
    regen_low,
    regen_high,
    regen_diff
  )
) %>%
  mutate(
    absolute_difference = abs(regenerated_value - manuscript_value),
    within_tolerance = absolute_difference <= MANUSCRIPT_CHECK_TOL_BPM
  )

write_csv(
  check_tbl,
  file.path(out_dir, "figure7_manuscript_reproducibility_check.csv")
)

log_msg(
  "Manuscript check | low=", sprintf("%.3f", regen_low),
  " | high=", sprintf("%.3f", regen_high),
  " | difference=", sprintf("%.3f", regen_diff)
)

if (!all(check_tbl$within_tolerance)) {
  warning(
    "Figure 7A/Table 5 values differ from the revised manuscript by more than ",
    MANUSCRIPT_CHECK_TOL_BPM,
    " bpm. Inspect figure7_manuscript_reproducibility_check.csv."
  )
} else {
  log_msg(
    "Figure 7A/Table 5 values agree with manuscript references within ",
    MANUSCRIPT_CHECK_TOL_BPM,
    " bpm."
  )
}

# ============================================================
# FIGURE 7A
# ============================================================

plot_a_df <- scenario_summary %>%
  mutate(
    scenario = factor(
      scenario,
      levels = names(SCENARIOS_HOURS_PER_DAY),
      labels = c(
        "30 min/day\ncommute",
        "2 hr/day\ncommute",
        "8 hr/day\nprofessional driving"
      )
    ),
    condition = factor(
      condition,
      levels = c(TRAIT_LABEL_LOW, TRAIT_LABEL_HIGH)
    )
  )

pA <- ggplot(
  plot_a_df,
  aes(x = scenario, y = annual_mean, fill = condition)
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.70,
    color = "black",
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(ymin = annual_q05, ymax = annual_q95),
    position = position_dodge(width = 0.78),
    width = 0.16,
    linewidth = 0.65
  ) +
  scale_fill_manual(values = PAL_TRAIT) +
  scale_y_continuous(
    labels = label_number(big.mark = ",", accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Annual cumulative NHR-hours\n[bpm·hours]",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(lineheight = 0.92),
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
    plot.margin = margin(6, 8, 4, 4)
  )

# ============================================================
# FIGURE 7B: HORIZON-SCALING FROM HELD-OUT ENET PREDICTIONS
# ============================================================

# Reproduce the former Figure 8 analysis exactly: reconstruct the row map
# from the clean MASTER dataset, retain common subjects across strata,
# join OOF predictions by stratum and row_id, build contiguous segments,
# and retain complete non-overlapping blocks only.
USE_COMMON_SUBJECTS <- TRUE
MIN_BLOCKS_WARN <- 30L

dt_full <- fread(data_path, showProgress = TRUE)
dt_full <- canonicalize_names(dt_full)

need_b <- c("p_id", "time", "raw_hr", "activity3", "bl_hr")
miss_b <- setdiff(need_b, names(dt_full))
if (length(miss_b) > 0L) {
  stop("MASTER data missing for Figure 7B: ", paste(miss_b, collapse = ", "))
}

dt_full[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]
dt_full <- dt_full[!is.na(dt_time)]
dt_full[, activity3_norm := norm_chr(activity3)]
dt_full[, day_key := make_day_key(dt_full)]

dt_drive_b <- dt_full[activity3_norm == "driving"]
dt_nond_b <- dt_full[activity3_norm == "non_driving_sedentary"]
if (nrow(dt_drive_b) == 0L) stop("No DRIVING rows found for Figure 7B.")
if (nrow(dt_nond_b) == 0L) stop("No NONDRIVING_SEDENTARY rows found for Figure 7B.")

if (USE_COMMON_SUBJECTS) {
  common_subj <- intersect(unique(dt_drive_b$p_id), unique(dt_nond_b$p_id))
  dt_drive_b <- dt_drive_b[p_id %in% common_subj]
  dt_nond_b <- dt_nond_b[p_id %in% common_subj]
  log_msg("Figure 7B common subjects: ", length(common_subj))
}

setorderv(dt_drive_b, c("p_id", "dt_time"))
setorderv(dt_nond_b, c("p_id", "dt_time"))

rowmap_drive <- copy(dt_drive_b)[
  , .(p_id, dt_time, day_key, stratum = "DRIVING")
][, row_id := seq_len(.N)]

rowmap_nond <- copy(dt_nond_b)[
  , .(p_id, dt_time, day_key, stratum = "NONDRIVING_SEDENTARY")
][, row_id := seq_len(.N)]

rowmap <- rbindlist(list(rowmap_drive, rowmap_nond), use.names = TRUE)
fwrite(rowmap_drive, file.path(out_dir, "figure7B_diagnostics_rowmap_DRIVING.csv"))
fwrite(rowmap_nond, file.path(out_dir, "figure7B_diagnostics_rowmap_NONDRIVING_SEDENTARY.csv"))

pred <- fread(run_inputs["predictions"])
need_pred <- c("model", "stratum", "fold", "row_id", "raw_hr_obs", "raw_hr_hat")
miss_pred <- setdiff(need_pred, names(pred))
if (length(miss_pred) > 0L) {
  stop("Prediction file missing columns: ", paste(miss_pred, collapse = ", "))
}

pred <- pred[model == "enet"]
if (nrow(pred) == 0L) stop("No ENet rows in prediction file.")

# The current prediction export may already contain participant/time mapping
# columns. Remove them before joining the authoritative row map reconstructed
# from the MASTER dataset; otherwise data.table::merge() creates p_id.x/p_id.y
# and dt_time.x/dt_time.y, leaving no column named p_id or dt_time.
mapping_cols_in_pred <- intersect(
  c("p_id", "dt_time", "day_key"),
  names(pred)
)
if (length(mapping_cols_in_pred) > 0L) {
  log_msg(
    "Dropping pre-existing prediction mapping columns before Figure 7B join: ",
    paste(mapping_cols_in_pred, collapse = ", ")
  )
  pred[, (mapping_cols_in_pred) := NULL]
}

pred <- merge(
  pred,
  rowmap,
  by = c("stratum", "row_id"),
  all.x = TRUE,
  all.y = FALSE,
  sort = FALSE
)

if (pred[, any(is.na(p_id) | is.na(dt_time))]) {
  bad_n <- pred[is.na(p_id) | is.na(dt_time), .N]
  stop("Failed to reconstruct row mapping for ", bad_n, " prediction rows.")
}

pred[, raw_resid := as.numeric(raw_hr_obs) - as.numeric(raw_hr_hat)]
setorderv(pred, c("stratum", "p_id", "dt_time"))

pred[, dt_diff_sec := as.numeric(
  difftime(dt_time, shift(dt_time), units = "secs")
), by = .(stratum, p_id)]
pred[, day_key_prev := shift(day_key), by = .(stratum, p_id)]

pred[, new_segment := fifelse(
  is.na(dt_diff_sec) |
    dt_diff_sec <= 0 |
    dt_diff_sec > MAX_GAP_MULTIPLIER * RES_SECONDS |
    is.na(day_key_prev) |
    day_key != day_key_prev,
  1L, 0L
), by = .(stratum, p_id)]

pred[, segment_id := cumsum(new_segment), by = .(stratum, p_id)]

compute_horizon_blocks <- function(dt_in, horizon_min, res_seconds) {
  k <- as.integer(round((horizon_min * 60) / res_seconds))
  if (k < 1L) stop("Invalid horizon.")

  copy(dt_in)[
    order(stratum, p_id, dt_time),
    {
      idx <- seq_len(.N)
      block_id <- ((idx - 1L) %/% k) + 1L
      .(
        horizon_min = horizon_min,
        block_id = block_id,
        raw_resid = raw_resid
      )
    },
    by = .(stratum, p_id, segment_id)
  ][
    ,
    .(
      n_rows = .N,
      raw_mean_resid = mean(raw_resid, na.rm = TRUE)
    ),
    by = .(stratum, p_id, segment_id, horizon_min, block_id)
  ][n_rows == k]
}

blocks_all <- rbindlist(
  lapply(HORIZONS_MIN, function(hm) {
    compute_horizon_blocks(pred, hm, RES_SECONDS)
  }),
  use.names = TRUE,
  fill = TRUE
)
if (nrow(blocks_all) == 0L) stop("No complete Figure 7B blocks formed.")

fwrite(blocks_all, file.path(out_dir, "figure7_horizon_block_means_long.csv"))

horizon_summary <- blocks_all[
  ,
  .(
    rmse = sqrt(mean(raw_mean_resid^2, na.rm = TRUE)),
    n_blocks = .N
  ),
  by = .(stratum, horizon_min)
]
horizon_summary[, normalized_rmse := rmse / rmse[horizon_min == min(horizon_min)], by = stratum]
horizon_summary[, warn_sparse := n_blocks < MIN_BLOCKS_WARN]
horizon_summary[, stratum_plot := factor(
  stratum,
  levels = c("DRIVING", "NONDRIVING_SEDENTARY"),
  labels = c("Driving", "Non-driving sedentary")
)]

write_csv(
  as_tibble(horizon_summary),
  file.path(out_dir, "figure7_horizon_scaling_summary.csv")
)

pB <- ggplot(
  horizon_summary,
  aes(
    x = horizon_min,
    y = normalized_rmse,
    color = stratum_plot,
    group = stratum_plot
  )
) +
  geom_line(linewidth = 1.25, lineend = "round") +
  geom_point(aes(shape = warn_sparse), size = 3.1, stroke = 1.15) +
  scale_color_manual(values = PAL_STRATUM, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_x_log10(
    breaks = HORIZONS_MIN,
    labels = HORIZONS_MIN,
    expand = expansion(mult = c(0.04, 0.04))
  ) +
  scale_y_log10(
    breaks = c(0.80, 0.90, 1.00),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  labs(
    x = "Averaging horizon [min, log scale]",
    y = "Normalized RMSE"
  ) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(0.85, "cm"),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(margin = margin(t = 9)),
    axis.title.y = element_text(margin = margin(r = 9)),
    axis.text = element_text(size = 12),
    axis.line = element_line(linewidth = 0.55, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    panel.grid.major = element_line(linewidth = 0.30, color = "grey88"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 4, 8)
  )

# ============================================================
# SAVE FIGURE 7
# ============================================================

fig7 <- (pA | pB) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

pdf_path <- file.path(fig_dir, "Figure7_LongHorizon_Behavior.pdf")
png_path <- file.path(fig_dir, "Figure7_LongHorizon_Behavior.png")

safe_save_pdf(fig7, pdf_path, w = 11.2, h = 5.4)
ggsave(
  png_path,
  fig7,
  width = 11.2,
  height = 5.4,
  dpi = PNG_DPI
)

log_msg("Saved: ", pdf_path)
log_msg("Saved: ", png_path)
log_msg("Wrote Table 5 support file and manuscript reproducibility check.")
log_msg("DONE. Outputs in: ", out_dir)
