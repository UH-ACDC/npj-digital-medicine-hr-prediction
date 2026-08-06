# ============================================================
# 05_figure5_table4_predictability_decomposition.R
#
# PURPOSE
#   Generate Figure 5 and the supporting Table 4 summary for the
#   npj Digital Public Health manuscript:
#
#     Wearable sensing reveals the structure of cardiac
#     activation associated with everyday driving
#
#   The script summarizes cross-validated prediction performance
#   for heart-rate prediction models in driving and non-driving
#   sedentary contexts.
#
# ANALYTIC DEFINITION
#   Figure 5 decomposes predictability into three nested model
#   stages:
#
#     baseline0
#       Person baseline only
#
#     baseline_offset
#       Person baseline plus learned context offset
#       - driving: driving-associated context offset
#       - non-driving sedentary: non-driving sedentary context offset
#
#     enet
#       Elastic-net model including additional modulators
#
# PANEL DEFINITIONS
#   a. Cross-validated RMSE across folds
#   b. Cross-validated squared Pearson correlation r^2 across folds
#
# INPUT
#   The script searches recursively under:
#
#     Results/
#
#   Preferred prediction filenames:
#
#     predictions_all_models_both_strata.csv
#     predictions_all_models.csv
#     oof_predictions_all_models_both_strata.csv
#     oof_predictions_all_models.csv
#
#   The script prioritizes files matching the selected resolution.
#
# MAJOR OUTPUTS
#   The script writes Figure 5, the Table 4 summary, and supporting
#   diagnostics under:
#
#     Results/paper_figs/<timestamp>_<RES>sec_figure5_table4_predictability_decomposition/
#
#   Main outputs:
#
#     Figures/Figure5_Predictability_Decomposition.pdf
#     Figures/Figure5_Predictability_Decomposition.png
#     figure5_metrics_by_fold.csv
#     figure5_summary.csv
#     figure5_decomposition_rmse.csv
#     table4_predictive_decomposition.csv
#     run_log.txt
#
# REPOSITORY SCOPE
#   This public repository starts from curated model prediction
#   outputs. It does not retrain the machine-learning models from
#   raw wearable, smartphone, vehicle, or ground-truth streams.
#
# PRIVACY NOTE
#   This script operates on fold-level prediction outputs and does
#   not require direct GPS coordinates or raw location traces.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggpattern)
  library(patchwork)
})

options(warn = 1)
set.seed(20260309)

# ----------------------------
# User toggles
# ----------------------------
PDF_W <- 10.5
PDF_H <- 8.0
PNG_DPI <- 300
SAVE_R2_IF_ALL_NA <- TRUE

MODEL_ORDER <- c("baseline0", "baseline_offset", "enet")

# Within-facet shades:
#   DRIVING = oranges
#   NONDRIVING_SEDENTARY = grays
FACET_MODEL_PAL <- c(
  "DRIVING__Baseline only"                             = "#FDD9B5",
  "DRIVING__Baseline + context offset"                = "#F4A261",
  "DRIVING__ENet"                                     = "#D97706",
  "NONDRIVING_SEDENTARY__Baseline only"               = "white",
  "NONDRIVING_SEDENTARY__Baseline + context offset"   = "grey80",
  "NONDRIVING_SEDENTARY__ENet"                        = "grey55"
)

PRED_BASENAMES <- c(
  "predictions_all_models_both_strata.csv",
  "predictions_all_models.csv",
  "oof_predictions_all_models_both_strata.csv",
  "oof_predictions_all_models.csv"
)

# ----------------------------
# Resolve paths RELATIVE to Scripts/
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
results_root <- normalizePath(file.path(project_root, "Results"), mustWork = TRUE)
paper_figs_root <- normalizePath(file.path(results_root, "paper_figs"), mustWork = FALSE)

if (!dir.exists(paper_figs_root)) {
  dir.create(paper_figs_root, recursive = TRUE, showWarnings = FALSE)
}

message("Script dir: ", script_dir)
message("Project root: ", project_root)
message("Results root: ", results_root)
message("paper_figs root: ", paper_figs_root)

# ----------------------------
# Analysis resolution
# ----------------------------
# The manuscript and public-repository outputs use 60-s data.
# Change this value manually only if matching 10-s or 30-s
# prediction files are added under Results/ in the future.
RES_SECONDS <- 60L

# ----------------------------
# Helpers
# ----------------------------
first_existing_col <- function(df, candidates, required = TRUE, what = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    if (required) {
      stop("Could not find ", what, ". Tried: ", paste(candidates, collapse = ", "))
    } else {
      return(NA_character_)
    }
  }
  hit[1]
}

normalize_stratum <- function(x) {
  x0 <- toupper(trimws(as.character(x)))
  x0 <- str_replace_all(x0, "[[:space:]/-]+", "_")
  
  dplyr::case_when(
    x0 %in% c("DRIVING", "DRIVE") ~ "DRIVING",
    x0 %in% c(
      "NONDRIVING_SEDENTARY",
      "NON_DRIVING_SEDENTARY",
      "NONDRIVING",
      "NON_DRIVING",
      "SEDENTARY_NONDRIVING",
      "SEDENTARY_NON_DRIVING",
      "NONDRIVINGSEDENTARY"
    ) ~ "NONDRIVING_SEDENTARY",
    TRUE ~ x0
  )
}

normalize_model <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- str_replace_all(x0, "[[:space:]-]+", "_")
  
  dplyr::case_when(
    x0 %in% c(
      "baseline0", "baseline_0", "baseline", "baseline_only",
      "baselineonly", "person_baseline", "participant_baseline"
    ) ~ "baseline0",
    x0 %in% c("baseline_offset", "baseline_plus_offset", "offset", "baselineoffset") ~ "baseline_offset",
    x0 %in% c("enet", "elastic_net", "elasticnet", "glmnet") ~ "enet",
    TRUE ~ x0
  )
}

extract_version_num <- function(x) {
  x_low <- tolower(x)
  m <- stringr::str_match(x_low, "v([0-9]+)")
  out <- suppressWarnings(as.numeric(m[, 2]))
  ifelse(is.na(out), -Inf, out)
}

rmse_vec <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) == 0L) return(NA_real_)
  sqrt(mean((y - p)^2))
}

rsq_vec <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) < 2L) return(NA_real_)
  if (sd(y) == 0 || sd(p) == 0) return(NA_real_)
  cor(y, p)^2
}

r2_oos_vec <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) < 2L) return(NA_real_)
  sst <- sum((y - mean(y))^2)
  if (!is.finite(sst) || sst <= 0) return(NA_real_)
  1 - sum((y - p)^2) / sst
}

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  sd(x) / sqrt(length(x))
}

safe_save_pdf <- function(plot_obj, path, w = 7, h = 5) {
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

safe_save_png <- function(plot_obj, path, w = 7, h = 5, dpi = 300) {
  ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}

# Stratum-specific display labels
model_label_by_stratum <- function(stratum, model) {
  dplyr::case_when(
    model == "baseline0" ~ "Baseline only",
    model == "baseline_offset" ~ "Baseline + context offset",
    model == "enet" ~ "ENet",
    TRUE ~ NA_character_
  )
}

# ----------------------------
# Find prediction file under Results/
# relative to Scripts/
# ----------------------------
find_prediction_candidates <- function(results_root, res_seconds) {
  all_csv <- list.files(
    results_root,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_csv) == 0) return(character(0))
  
  base_ok <- basename(all_csv) %in% PRED_BASENAMES
  
  res_pat <- paste0("(^|[^0-9])", res_seconds, "sec([^0-9]|$)")
  res_ok <- grepl(res_pat, all_csv, ignore.case = TRUE, perl = TRUE)
  
  cand <- all_csv[base_ok & res_ok]
  if (length(cand) == 0) cand <- all_csv[base_ok]
  
  cand
}

prediction_candidates <- find_prediction_candidates(results_root, RES_SECONDS)

if (length(prediction_candidates) == 0) {
  stop(
    "Could not find any prediction CSV under: ", results_root, "\n",
    "Looked for basenames: ", paste(PRED_BASENAMES, collapse = ", ")
  )
}

candidate_rank <- data.frame(
  path = prediction_candidates,
  is_both = grepl("both_strata", basename(prediction_candidates), ignore.case = TRUE),
  version_num = extract_version_num(prediction_candidates),
  mtime = file.info(prediction_candidates)$mtime,
  stringsAsFactors = FALSE
) %>%
  arrange(desc(version_num), desc(is_both), desc(mtime))

pred_path <- candidate_rank$path[1]
run_dir <- dirname(pred_path)

message("Selected prediction file: ", pred_path)
message("Inferred run dir: ", run_dir)

# ----------------------------
# Output folder under Results/paper_figs
# ----------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_dir <- file.path(
  paper_figs_root,
  paste0(stamp, "_", RES_SECONDS, "sec_figure5_table4_predictability_decomposition")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("Script: 05_figure5_table4_predictability_decomposition.R")
log_msg("Figure 5 / Table 4 predictive decomposition analysis start")
log_msg("Project root: ", project_root)
log_msg("Results root: ", results_root)
log_msg("Resolution: ", RES_SECONDS, " sec")
log_msg("Run dir inferred: ", run_dir)
log_msg("Pred path: ", pred_path)
log_msg("Output: ", out_dir)
log_msg("Candidate ranking:")
capture.output(print(candidate_rank), file = log_file, append = TRUE)

# ============================================================
# Load predictions
# ============================================================
preds_raw <- suppressWarnings(readr::read_csv(pred_path, show_col_types = FALSE))
log_msg("Loaded predictions with n = ", nrow(preds_raw), " rows and p = ", ncol(preds_raw), " columns")

col_model   <- first_existing_col(preds_raw, c("model", ".model", "model_name"), TRUE, "model column")
col_stratum <- first_existing_col(preds_raw, c("stratum", "context", "subset", "activity3", "activity"), TRUE, "stratum column")
col_fold    <- first_existing_col(preds_raw, c("fold", "id", "resample", "cv_fold", "fold_id"), TRUE, "fold column")
col_obs     <- first_existing_col(preds_raw, c("raw_hr_obs", "obs", "truth", "y", ".outcome"), TRUE, "observed outcome column")
col_hat     <- first_existing_col(preds_raw, c("raw_hr_hat", "pred", ".pred", "estimate", "prediction", "yhat"), TRUE, "prediction column")

log_msg("Raw model labels: ", paste(sort(unique(as.character(preds_raw[[col_model]]))), collapse = ", "))

log_msg("Detected model col: ", col_model)
log_msg("Detected stratum col: ", col_stratum)
log_msg("Detected fold col: ", col_fold)
log_msg("Detected obs col: ", col_obs)
log_msg("Detected pred col: ", col_hat)

preds <- preds_raw %>%
  transmute(
    model      = normalize_model(.data[[col_model]]),
    stratum    = normalize_stratum(.data[[col_stratum]]),
    fold       = as.character(.data[[col_fold]]),
    raw_hr_obs = suppressWarnings(as.numeric(.data[[col_obs]])),
    raw_hr_hat = suppressWarnings(as.numeric(.data[[col_hat]]))
  ) %>%
  filter(model %in% MODEL_ORDER) %>%
  filter(stratum %in% c("DRIVING", "NONDRIVING_SEDENTARY"))

if (nrow(preds) == 0) {
  stop("No rows remain after harmonizing models/strata. Check your predictions file: ", pred_path)
}

log_msg("Rows after filtering to paper models/strata: ", nrow(preds))
log_msg("Models present: ", paste(sort(unique(preds$model)), collapse = ", "))
log_msg("Strata present: ", paste(sort(unique(preds$stratum)), collapse = ", "))

missing_models <- setdiff(MODEL_ORDER, unique(preds$model))
if (length(missing_models) > 0L) {
  stop(
    "Required model stage(s) missing after label harmonization: ",
    paste(missing_models, collapse = ", "), "\n",
    "Raw model labels in the selected prediction file: ",
    paste(sort(unique(as.character(preds_raw[[col_model]]))), collapse = ", ")
  )
}

# ============================================================
# Fold-level metrics
# ============================================================
by_fold <- preds %>%
  group_by(stratum, model, fold) %>%
  summarise(
    n    = sum(is.finite(raw_hr_obs) & is.finite(raw_hr_hat)),
    rmse   = rmse_vec(raw_hr_obs, raw_hr_hat),
    rsq    = rsq_vec(raw_hr_obs, raw_hr_hat),
    r2_oos = r2_oos_vec(raw_hr_obs, raw_hr_hat),
    .groups = "drop"
  ) %>%
  arrange(stratum, match(model, MODEL_ORDER), fold)

write_csv(by_fold, file.path(out_dir, "figure5_metrics_by_fold.csv"))
log_msg("Wrote figure5_metrics_by_fold.csv")

# ============================================================
# Summary across folds
# ============================================================
sum_fold <- by_fold %>%
  group_by(stratum, model) %>%
  summarise(
    k_folds   = sum(is.finite(rmse)),
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_se   = safe_se(rmse),
    rsq_mean    = mean(rsq, na.rm = TRUE),
    rsq_se      = safe_se(rsq),
    r2_oos_mean = mean(r2_oos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    stratum = factor(stratum, levels = c("DRIVING", "NONDRIVING_SEDENTARY")),
    model   = factor(model, levels = MODEL_ORDER)
  ) %>%
  arrange(stratum, model)

write_csv(sum_fold, file.path(out_dir, "figure5_summary.csv"))
log_msg("Wrote figure5_summary.csv")

# ============================================================
# Decomposition table for RMSE narration
# ============================================================
decomp_rmse <- sum_fold %>%
  mutate(
    stratum = as.character(stratum),
    model   = as.character(model)
  ) %>%
  select(stratum, model, rmse_mean) %>%
  pivot_wider(names_from = model, values_from = rmse_mean) %>%
  mutate(
    improv_offset     = baseline0 - baseline_offset,
    improv_modulators = baseline_offset - enet,
    total_improv      = baseline0 - enet,
    frac_offset = ifelse(
      is.finite(total_improv) & total_improv != 0,
      improv_offset / total_improv,
      NA_real_
    ),
    frac_modulators = ifelse(
      is.finite(total_improv) & total_improv != 0,
      improv_modulators / total_improv,
      NA_real_
    )
  )

write_csv(decomp_rmse, file.path(out_dir, "figure5_decomposition_rmse.csv"))
log_msg("Wrote figure5_decomposition_rmse.csv")

# ============================================================
# Table 4 summary
# ============================================================
# Delta metrics are computed relative to the baseline-only model.
# "ENet wins" is the percentage of outer folds in which ENet has
# lower RMSE than the baseline-plus-offset model.

table4 <- sum_fold %>%
  mutate(
    stratum = as.character(stratum),
    model = as.character(model)
  ) %>%
  group_by(stratum) %>%
  mutate(
    delta_rmse = rmse_mean - rmse_mean[model == "baseline0"],
    delta_rsq  = rsq_mean - rsq_mean[model == "baseline0"]
  ) %>%
  ungroup()

enet_wins <- by_fold %>%
  select(stratum, model, fold, rmse) %>%
  filter(model %in% c("baseline_offset", "enet")) %>%
  pivot_wider(names_from = model, values_from = rmse) %>%
  group_by(stratum) %>%
  summarise(
    enet_wins_pct = 100 * mean(enet < baseline_offset, na.rm = TRUE),
    .groups = "drop"
  )

table4 <- table4 %>%
  left_join(enet_wins, by = "stratum") %>%
  mutate(
    model_label = case_when(
      model == "baseline0" ~ "Baseline only",
      model == "baseline_offset" ~ "Baseline + context offset",
      model == "enet" ~ "ENet",
      TRUE ~ model
    ),
    r2_oos_report = ifelse(model == "enet", r2_oos_mean, NA_real_),
    enet_wins_report = ifelse(model == "enet", enet_wins_pct, NA_real_)
  ) %>%
  select(
    stratum,
    model,
    model_label,
    k_folds,
    rmse_mean,
    rmse_se,
    delta_rmse,
    rsq_mean,
    rsq_se,
    delta_rsq,
    r2_oos_report,
    enet_wins_report
  ) %>%
  arrange(match(stratum, c("DRIVING", "NONDRIVING_SEDENTARY")), match(model, MODEL_ORDER))

write_csv(table4, file.path(out_dir, "table4_predictive_decomposition.csv"))
log_msg("Wrote table4_predictive_decomposition.csv")

# ============================================================
# Plot prep
# ============================================================
plot_df <- sum_fold %>%
  mutate(
    stratum_raw = as.character(stratum),
    stratum = factor(
      stratum_raw,
      levels = c("DRIVING", "NONDRIVING_SEDENTARY"),
      labels = c("Driving", "Non-driving sedentary")
    ),
    model_lab = model_label_by_stratum(stratum_raw, as.character(model))
  ) %>%
  mutate(
    model_lab = factor(
      model_lab,
      levels = c(
        "Baseline only",
        "Baseline + context offset",
        "ENet"
      )
    ),
    palette_key = paste(stratum_raw, as.character(model_lab), sep = "__"),
    fill_col = unname(FACET_MODEL_PAL[palette_key]),
    fill_col = ifelse(is.na(fill_col), "white", fill_col),
    pattern_type = case_when(
      stratum_raw == "NONDRIVING_SEDENTARY" &
        model_lab == "Baseline only" ~ "stripe",
      TRUE ~ "none"
    ),
    err_col = "black"
  )

theme_fig5 <- theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11, angle = 12, hjust = 1),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11)
  )

# ============================================================
# Figure 5a: RMSE
# ============================================================
p_rmse <- ggplot(plot_df, aes(x = model_lab, y = rmse_mean)) +
  geom_col_pattern(
    aes(
      fill = fill_col,
      pattern = pattern_type
    ),
    width = 0.72,
    color = "black",
    linewidth = 0.25,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.08,
    pattern_spacing = 0.03,
    pattern_size = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = rmse_mean - rmse_se,
      ymax = rmse_mean + rmse_se,
      color = err_col
    ),
    width = 0.18,
    linewidth = 0.8,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_pattern_identity() +
  scale_color_identity() +
  facet_wrap(~ stratum, nrow = 1, scales = "free_x") +
  labs(
    subtitle = "Cross-validated RMSE across outer folds",
    x = NULL,
    y = "RMSE [bpm]"
  ) +
  theme_fig5

# ============================================================
# Figure 5b: R^2
# ============================================================
p_rsq <- ggplot(plot_df, aes(x = model_lab, y = rsq_mean)) +
  geom_col_pattern(
    aes(
      fill = fill_col,
      pattern = pattern_type
    ),
    width = 0.72,
    color = "black",
    linewidth = 0.25,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.08,
    pattern_spacing = 0.03,
    pattern_size = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = rsq_mean - rsq_se,
      ymax = rsq_mean + rsq_se,
      color = err_col
    ),
    width = 0.18,
    linewidth = 0.8,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_pattern_identity() +
  scale_color_identity() +
  facet_wrap(~ stratum, nrow = 1, scales = "free_x") +
  labs(
    subtitle = "Cross-validated variance explained across outer folds",
    x = NULL,
    y = expression(r^2)
  ) +
  theme_fig5

# ============================================================
# Assemble and save composite Figure 5
# ============================================================
# Springer Nature requests all panels belonging to one figure to be
# supplied in a single file. Lowercase panel tags are used to match
# manuscript references such as Fig. 5a and Fig. 5b.

all_rsq_na <- all(!is.finite(plot_df$rsq_mean))

if (all_rsq_na && !isTRUE(SAVE_R2_IF_ALL_NA)) {
  stop(
    "Figure 5b cannot be assembled because all R2 values are NA and ",
    "SAVE_R2_IF_ALL_NA is FALSE."
  )
}

# Add a little more whitespace between the upper and lower panels.
p_rmse <- p_rmse +
  theme(plot.margin = margin(t = 5, r = 5, b = 14, l = 5))

p_rsq <- p_rsq +
  theme(plot.margin = margin(t = 14, r = 5, b = 5, l = 5))

figure5 <- (p_rmse / p_rsq) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.01, 0.99)
  )

safe_save_pdf(
  figure5,
  file.path(fig_dir, "Figure5_Predictability_Decomposition.pdf"),
  w = PDF_W,
  h = PDF_H
)

safe_save_png(
  figure5,
  file.path(fig_dir, "Figure5_Predictability_Decomposition.png"),
  w = PDF_W,
  h = PDF_H,
  dpi = PNG_DPI
)

log_msg("Saved composite Figure 5 with panels a and b")
log_msg("Saved figures to: ", fig_dir)

# ============================================================
# Log headline decomposition numbers
# ============================================================
log_msg("Decomposition (RMSE drops):")
capture.output(print(decomp_rmse), file = log_file, append = TRUE)

log_msg("DONE. Figure 5 and Table 4 outputs in: ", out_dir)