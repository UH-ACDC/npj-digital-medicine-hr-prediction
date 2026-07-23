# ============================================================
# 08_supplementary_tableS1_hr_missingness_analysis_v1.1.R
#
# PURPOSE
#   Reproduce the heart-rate missingness analysis reported in the Supplementary
#   Materials and summarized in Supplementary Table S1. The analysis quantifies
#   whether raw heart-rate measurements are more frequently missing during
#   driving than during non-driving sedentary behavior and, within driving,
#   whether missingness is associated with objective indicators of driving demand.
#
# INPUT
#   Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
# OUTPUTS
#   Results are written to:
#     Results/paper_figs/<timestamp>_60sec_supplementary_tableS1_hr_missingness/
#
#   Core Supplement/Table S1 outputs:
#     missingness_driving_vs_sedentary_summary.csv
#     missingness_model_primary_driving_vs_sedentary_odds_ratios.csv
#     missingness_model_driving_intensity_focused_odds_ratios.csv
#
#   Supporting descriptive, sensitivity, audit, and diagnostic outputs are also
#   written to the same directory. See missingness_analysis_notes.txt for the
#   complete model specifications and output map.
#
# ANALYTIC STRUCTURE
#   1) Supplementary primary model (reported in Table S1):
#        HR_missing ~ activity_binary + day_period + day_type + weather_info
#                     + (1 | p_id)
#      Estimates the adjusted association between behavioral context and raw-HR
#      missingness while accounting for repeated observations within participants.
#
#   2) All-activity sensitivity model:
#        HR_missing ~ activity3 + day_period + day_type + weather_info
#                     + participant traits + (1 | p_id)
#      This estimates whether driving has higher missingness than non-driving
#      sedentary behavior after adjusting for broad temporal/contextual factors.
#
#   3) Supplementary driving-only model (reported in Table S1):
#        HR_missing ~ speed + jf + energy_acc + energy_rot + (1 | p_id)
#      Tests whether missingness within driving is concentrated in epochs with
#      higher speed, congestion, or wrist-motion energy.
#
#   4) Adjusted driving-only sensitivity model:
#        HR_missing ~ weather_info + day_period + day_type + driving indicators
#                     + participant traits + (1 | p_id)
#      Evaluates the stability of the driving-only associations after broader
#      temporal, environmental, and participant-level adjustment.
#
# NOTES
#   - Mixed-effects logistic models are fit with lme4::glmer().
#   - Numeric predictors are z-scored; corresponding odds ratios represent a
#     one-standard-deviation increase.
#   - Additional descriptive and adjusted sensitivity outputs are written to the
#     same timestamped directory for auditability but are not required for Table S1.
# ============================================================

options(warn = 1)
set.seed(20260701)

# ----------------------------
# User settings
# ----------------------------
# Robust project-root detection. This allows the script to be run either from
# the project root, the Scripts/ folder, or an arbitrary working directory after
# setting PROJECT_ROOT manually.
PROJECT_ROOT <- Sys.getenv("NUBI_PROJECT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT) || PROJECT_ROOT == "") {
  candidate_roots <- unique(normalizePath(c(getwd(), file.path(getwd(), "..")), mustWork = FALSE))
  has_data_file <- file.exists(file.path(candidate_roots, "Data", "NUBI_Data_60sec_Level_MASTER_CLEAN.csv"))
  if (!any(has_data_file)) {
    stop(
      "Could not locate Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv. ",
      "Run this script from the project root or Scripts/ folder, or set ",
      "Sys.setenv(NUBI_PROJECT_ROOT = '/path/to/project')."
    )
  }
  PROJECT_ROOT <- candidate_roots[which(has_data_file)[1]]
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = TRUE)

DATA_FILE <- file.path(
  PROJECT_ROOT,
  "Data",
  "NUBI_Data_60sec_Level_MASTER_CLEAN.csv"
)

OUT_ROOT <- file.path(PROJECT_ROOT, "Results", "paper_figs")
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

MIN_WEATHER_N <- 100

# ----------------------------
# Packages
# ----------------------------
required <- c("dplyr", "readr", "ggplot2", "forcats", "stringr", "tibble", "tidyr", "scales")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them before running this script.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(scales)
})

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Missing required package: lme4\nInstall it before running this script.")
}
HAS_LME4 <- TRUE

# ----------------------------
# Output folders
# ----------------------------
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(OUT_ROOT, paste0(stamp, "_60sec_supplementary_tableS1_hr_missingness"))
FIG_DIR <- file.path(OUT_DIR, "Figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Helper functions
# ----------------------------
check_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - m) / s)
}

summarize_missingness <- function(data, ...) {
  data %>%
    group_by(...) %>%
    summarise(
      n_rows = n(),
      n_missing_hr = sum(hr_missing, na.rm = TRUE),
      missing_rate = mean(hr_missing, na.rm = TRUE),
      n_subjects = n_distinct(p_id),
      .groups = "drop"
    ) %>%
    mutate(
      missing_percent = 100 * missing_rate,
      ci_low = 100 * pmax(0, missing_rate - 1.96 * sqrt(missing_rate * (1 - missing_rate) / n_rows)),
      ci_high = 100 * pmin(1, missing_rate + 1.96 * sqrt(missing_rate * (1 - missing_rate) / n_rows))
    )
}

extract_or_table <- function(model, model_label, fallback = FALSE) {
  sm <- summary(model)
  coefs <- as.data.frame(sm$coefficients)
  coefs$term <- rownames(coefs)
  rownames(coefs) <- NULL

  # Exact extraction avoids the earlier bug where p_value could be copied from
  # the wrong column because of regex alternation involving | symbols.
  estimate_col <- intersect(c("Estimate"), names(coefs))[1]
  se_col <- intersect(c("Std. Error"), names(coefs))[1]
  stat_col <- intersect(c("z value", "t value"), names(coefs))[1]
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), names(coefs))[1]

  if (is.na(estimate_col) || is.na(se_col)) {
    stop("Could not identify Estimate and Std. Error columns in model summary.")
  }

  out <- coefs %>%
    transmute(
      model = model_label,
      term = term,
      estimate_log_odds = .data[[estimate_col]],
      std_error = .data[[se_col]],
      statistic = if (!is.na(stat_col)) .data[[stat_col]] else NA_real_,
      p_value = if (!is.na(p_col)) .data[[p_col]] else NA_real_,
      odds_ratio = exp(estimate_log_odds),
      conf_low = exp(estimate_log_odds - 1.96 * std_error),
      conf_high = exp(estimate_log_odds + 1.96 * std_error),
      fit_type = ifelse(fallback, "glm_with_participant_fixed_effects", "mixed_effects_logistic")
    ) %>%
    filter(term != "(Intercept)")

  out
}

fit_missingness_model <- function(formula_mixed, formula_fallback, data, model_label) {
  notes <- character()
  model <- NULL
  fallback <- FALSE

  if (HAS_LME4) {
    notes <- c(notes, "Attempted mixed-effects logistic regression using lme4::glmer().")
    model <- tryCatch(
      lme4::glmer(
        formula_mixed,
        data = data,
        family = stats::binomial(link = "logit"),
        nAGQ = 0,
        control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      ),
      error = function(e) {
        notes <<- c(notes, paste("glmer failed:", conditionMessage(e)))
        NULL
      },
      warning = function(w) {
        notes <<- c(notes, paste("glmer warning:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    )
  } else {
    notes <- c(notes, "lme4 not installed; skipped mixed-effects model.")
  }

  if (is.null(model)) {
    fallback <- TRUE
    notes <- c(notes, "Falling back to ordinary logistic regression with participant fixed effects.")
    model <- stats::glm(formula_fallback, data = data, family = stats::binomial(link = "logit"))
  }

  or_tab <- extract_or_table(model, model_label, fallback = fallback)
  list(model = model, or_table = or_tab, notes = notes, fallback = fallback)
}

clean_term_labels <- function(x) {
  x %>%
    str_replace_all("activity_binary", "Activity: ") %>%
    str_replace_all("activity3", "Activity: ") %>%
    str_replace_all("weather_info", "Weather: ") %>%
    str_replace_all("day_period", "Period: ") %>%
    str_replace_all("day_type", "Day type: ") %>%
    str_replace_all("gender", "Gender: ") %>%
    str_replace_all("speed_z", "Speed (per SD)") %>%
    str_replace_all("jf_z", "Jam factor (per SD)") %>%
    str_replace_all("energy_acc_z", "Acceleration energy (per SD)") %>%
    str_replace_all("energy_rot_z", "Rotational energy (per SD)") %>%
    str_replace_all("age_z", "Age (per SD)") %>%
    str_replace_all("trait_anxiety_z", "Trait anxiety (per SD)") %>%
    str_replace_all("extraversion_z", "Extraversion (per SD)") %>%
    str_replace_all("agreeableness_z", "Agreeableness (per SD)") %>%
    str_replace_all("conscientiousness_z", "Conscientiousness (per SD)") %>%
    str_replace_all("neuroticism_z", "Neuroticism (per SD)") %>%
    str_replace_all("openness_z", "Openness (per SD)")
}

# ----------------------------
# Load and prepare data
# ----------------------------
df <- readr::read_csv(DATA_FILE, show_col_types = FALSE)

check_cols(df, c("p_id", "activity3", "raw_hr", "day_period", "day_type", "weather_info"))

dat <- df %>%
  mutate(
    p_id = factor(p_id),
    hr_missing = as.integer(is.na(raw_hr)),
    activity3 = factor(activity3),
    activity3 = forcats::fct_relevel(activity3, "non_driving_sedentary"),
    day_period = factor(day_period),
    day_period = forcats::fct_relevel(day_period, "morning"),
    day_type = factor(day_type),
    day_type = forcats::fct_relevel(day_type, "weekdays"),
    weather_info = ifelse(is.na(weather_info) | weather_info == "", "missing", as.character(weather_info)),
    weather_info = forcats::fct_lump_min(factor(weather_info), min = MIN_WEATHER_N, other_level = "other"),
    activity_binary = dplyr::case_when(
      activity3 == "non_driving_sedentary" ~ "non_driving_sedentary",
      activity3 == "driving" ~ "driving",
      TRUE ~ NA_character_
    ),
    activity_binary = factor(activity_binary, levels = c("non_driving_sedentary", "driving"))
  )

# Relevel weather only if clear exists after lumping.
if ("clear" %in% levels(dat$weather_info)) {
  dat$weather_info <- forcats::fct_relevel(dat$weather_info, "clear")
}

# Add scaled numeric covariates when available.
for (v in c("speed", "ff_speed", "jf", "energy_acc", "energy_rot", "age", "trait_anxiety",
            "extraversion", "agreeableness", "conscientiousness", "neuroticism", "openness")) {
  if (v %in% names(dat)) dat[[paste0(v, "_z")]] <- safe_scale(dat[[v]])
}

if ("gender" %in% names(dat)) dat$gender <- factor(dat$gender)

# ----------------------------
# Descriptive tables
# ----------------------------
row_level_summary <- tibble(
  n_rows = nrow(dat),
  n_subjects = n_distinct(dat$p_id),
  n_missing_hr = sum(dat$hr_missing),
  missing_rate = mean(dat$hr_missing),
  missing_percent = 100 * missing_rate
)
write_csv(row_level_summary, file.path(OUT_DIR, "missingness_row_level_summary.csv"))

activity_summary <- summarize_missingness(dat, activity3) %>%
  arrange(desc(missing_percent))
write_csv(activity_summary, file.path(OUT_DIR, "missingness_activity_summary.csv"))

activity_time_summary <- summarize_missingness(dat, activity3, day_period, day_type) %>%
  arrange(activity3, day_type, day_period)
write_csv(activity_time_summary, file.path(OUT_DIR, "missingness_activity_by_dayperiod_daytype.csv"))

# Focused descriptive comparison for the supplementary two-context comparison.
activity_binary_summary <- dat %>%
  filter(!is.na(activity_binary)) %>%
  summarize_missingness(activity_binary) %>%
  arrange(activity_binary)
write_csv(activity_binary_summary, file.path(OUT_DIR, "missingness_driving_vs_sedentary_summary.csv"))

# Driving-only context bins for observed descriptive checks.
driving_dat <- dat %>% filter(activity3 == "driving")

make_quantile_bin <- function(x, label) {
  if (all(is.na(x))) return(factor(rep(NA_character_, length(x))))
  qs <- unique(stats::quantile(x, probs = seq(0, 1, 0.25), na.rm = TRUE, type = 7))
  if (length(qs) < 3) return(factor(rep("not enough variation", length(x))))
  cut(x, breaks = qs, include.lowest = TRUE, dig.lab = 4,
      labels = paste0(label, " Q", seq_len(length(qs) - 1)))
}

context_bins <- driving_dat %>%
  mutate(
    speed_bin = if ("speed" %in% names(.)) make_quantile_bin(speed, "Speed") else factor(NA_character_),
    jf_bin = if ("jf" %in% names(.)) make_quantile_bin(jf, "Jam factor") else factor(NA_character_),
    energy_acc_bin = if ("energy_acc" %in% names(.)) make_quantile_bin(energy_acc, "Accel energy") else factor(NA_character_),
    energy_rot_bin = if ("energy_rot" %in% names(.)) make_quantile_bin(energy_rot, "Rot energy") else factor(NA_character_)
  )

context_summary <- bind_rows(
  summarize_missingness(context_bins, weather_info) %>% mutate(context = "weather_info", level = as.character(weather_info)) %>% select(context, level, everything(), -weather_info),
  summarize_missingness(context_bins, speed_bin) %>% mutate(context = "speed_quartile", level = as.character(speed_bin)) %>% select(context, level, everything(), -speed_bin),
  summarize_missingness(context_bins, jf_bin) %>% mutate(context = "jam_factor_quartile", level = as.character(jf_bin)) %>% select(context, level, everything(), -jf_bin),
  summarize_missingness(context_bins, energy_acc_bin) %>% mutate(context = "acceleration_energy_quartile", level = as.character(energy_acc_bin)) %>% select(context, level, everything(), -energy_acc_bin),
  summarize_missingness(context_bins, energy_rot_bin) %>% mutate(context = "rotational_energy_quartile", level = as.character(energy_rot_bin)) %>% select(context, level, everything(), -energy_rot_bin)
) %>%
  filter(!is.na(level))
write_csv(context_summary, file.path(OUT_DIR, "missingness_driving_context_bins.csv"))

# ----------------------------
# Model 1: primary driving-vs-sedentary comparison
# ----------------------------
primary_covariates <- c("activity_binary", "day_period", "day_type", "weather_info")

primary_dat_pre <- dat %>%
  filter(!is.na(activity_binary)) %>%
  select(hr_missing, p_id, all_of(primary_covariates))

model_primary_dat <- primary_dat_pre %>% tidyr::drop_na()

form_primary_mixed <- as.formula(paste("hr_missing ~", paste(primary_covariates, collapse = " + "), "+ (1 | p_id)"))
form_primary_glm   <- as.formula(paste("hr_missing ~", paste(primary_covariates, collapse = " + "), "+ p_id"))

fit_primary <- fit_missingness_model(
  form_primary_mixed,
  form_primary_glm,
  model_primary_dat,
  "primary_driving_vs_sedentary"
)
write_csv(fit_primary$or_table, file.path(OUT_DIR, "missingness_model_primary_driving_vs_sedentary_odds_ratios.csv"))
writeLines(fit_primary$notes, file.path(OUT_DIR, "missingness_model_primary_fit_notes.txt"))

# ----------------------------
# Model 2: all epochs sensitivity model
# ----------------------------
all_covariates <- c("activity3", "day_period", "day_type", "weather_info")
for (v in c("age_z", "trait_anxiety_z", "extraversion_z", "agreeableness_z", "conscientiousness_z", "neuroticism_z", "openness_z", "gender")) {
  if (v %in% names(dat)) all_covariates <- c(all_covariates, v)
}

model_all_dat <- dat %>%
  select(hr_missing, p_id, all_of(all_covariates)) %>%
  tidyr::drop_na()

form_all_mixed <- as.formula(paste("hr_missing ~", paste(all_covariates, collapse = " + "), "+ (1 | p_id)"))
form_all_glm   <- as.formula(paste("hr_missing ~", paste(all_covariates, collapse = " + "), "+ p_id"))

fit_all <- fit_missingness_model(form_all_mixed, form_all_glm, model_all_dat, "all_epochs")
write_csv(fit_all$or_table, file.path(OUT_DIR, "missingness_model_all_epochs_odds_ratios.csv"))
writeLines(fit_all$notes, file.path(OUT_DIR, "missingness_model_all_epochs_fit_notes.txt"))

# ----------------------------
# Model 3: focused driving-intensity model
# ----------------------------
# This is the most direct model for the supplementary driving-only question: among driving
# epochs, is HR missingness associated with objective indicators of driving
# demand or steering-related wrist movement?
driving_intensity_covariates <- character(0)
for (v in c("speed_z", "jf_z", "energy_acc_z", "energy_rot_z")) {
  if (v %in% names(driving_dat)) driving_intensity_covariates <- c(driving_intensity_covariates, v)
}

if (length(driving_intensity_covariates) == 0) {
  stop("No driving-intensity covariates were found. Expected one or more of speed_z, jf_z, energy_acc_z, energy_rot_z.")
}

model_driving_intensity_dat <- driving_dat %>%
  select(hr_missing, p_id, all_of(driving_intensity_covariates)) %>%
  tidyr::drop_na()

form_driving_intensity_mixed <- as.formula(paste("hr_missing ~", paste(driving_intensity_covariates, collapse = " + "), "+ (1 | p_id)"))
form_driving_intensity_glm   <- as.formula(paste("hr_missing ~", paste(driving_intensity_covariates, collapse = " + "), "+ p_id"))

# ----------------------------
# Model 4: adjusted driving-only context sensitivity model
# ----------------------------
# This fuller model keeps weather, time context, and participant traits as
# adjustment covariates. It is useful as a sensitivity analysis, but the
# focused driving-intensity model above is the cleaner reported driving-only result.
driving_covariates <- c("weather_info", "day_period", "day_type")
for (v in c("speed_z", "jf_z", "energy_acc_z", "energy_rot_z", "age_z", "trait_anxiety_z",
            "extraversion_z", "agreeableness_z", "conscientiousness_z", "neuroticism_z", "openness_z", "gender")) {
  if (v %in% names(driving_dat)) driving_covariates <- c(driving_covariates, v)
}

model_driving_dat <- driving_dat %>%
  select(hr_missing, p_id, all_of(driving_covariates)) %>%
  tidyr::drop_na()

form_driving_mixed <- as.formula(paste("hr_missing ~", paste(driving_covariates, collapse = " + "), "+ (1 | p_id)"))
form_driving_glm   <- as.formula(paste("hr_missing ~", paste(driving_covariates, collapse = " + "), "+ p_id"))

# Audit how much data are retained after model-specific complete-case filtering.
retention_row <- function(label, before, after) {
  tibble(
    model = label,
    n_rows_before = nrow(before),
    n_rows_after = nrow(after),
    rows_retained_percent = 100 * n_rows_after / n_rows_before,
    n_subjects_before = n_distinct(before$p_id),
    n_subjects_after = n_distinct(after$p_id),
    missing_percent_before = 100 * mean(before$hr_missing, na.rm = TRUE),
    missing_percent_after = 100 * mean(after$hr_missing, na.rm = TRUE)
  )
}
row_retention_audit <- bind_rows(
  retention_row("primary_driving_vs_sedentary", primary_dat_pre, model_primary_dat),
  retention_row("all_epochs_sensitivity", dat %>% select(hr_missing, p_id, all_of(all_covariates)), model_all_dat),
  retention_row("driving_only_intensity_focused", driving_dat %>% select(hr_missing, p_id, all_of(driving_intensity_covariates)), model_driving_intensity_dat),
  retention_row("driving_only_context_adjusted", driving_dat %>% select(hr_missing, p_id, all_of(driving_covariates)), model_driving_dat)
)
write_csv(row_retention_audit, file.path(OUT_DIR, "missingness_model_row_retention_audit.csv"))

fit_driving_intensity <- fit_missingness_model(
  form_driving_intensity_mixed,
  form_driving_intensity_glm,
  model_driving_intensity_dat,
  "driving_only_intensity_focused"
)
write_csv(fit_driving_intensity$or_table, file.path(OUT_DIR, "missingness_model_driving_intensity_focused_odds_ratios.csv"))
writeLines(fit_driving_intensity$notes, file.path(OUT_DIR, "missingness_model_driving_intensity_focused_fit_notes.txt"))

fit_driving <- fit_missingness_model(form_driving_mixed, form_driving_glm, model_driving_dat, "driving_only_context_adjusted")
write_csv(fit_driving$or_table, file.path(OUT_DIR, "missingness_model_driving_context_adjusted_odds_ratios.csv"))
writeLines(fit_driving$notes, file.path(OUT_DIR, "missingness_model_driving_context_adjusted_fit_notes.txt"))

# ----------------------------
# Predicted probabilities for primary driving-vs-sedentary model
# ----------------------------
new_primary <- model_primary_dat %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  slice(rep(1, length(levels(model_primary_dat$activity_binary))))
new_primary$activity_binary <- factor(levels(model_primary_dat$activity_binary), levels = levels(model_primary_dat$activity_binary))
new_primary$p_id <- model_primary_dat$p_id[1]
for (v in c("day_period", "day_type", "weather_info")) {
  if (v %in% names(model_primary_dat)) {
    ref <- levels(model_primary_dat[[v]])[1]
    new_primary[[v]] <- factor(ref, levels = levels(model_primary_dat[[v]]))
  }
}

pred_primary <- tryCatch({
  if (!fit_primary$fallback && HAS_LME4) {
    p <- predict(fit_primary$model, newdata = new_primary, type = "response", re.form = NA, allow.new.levels = TRUE)
  } else {
    p <- predict(fit_primary$model, newdata = new_primary, type = "response")
  }
  tibble(activity_binary = new_primary$activity_binary, predicted_missing_rate = p, predicted_missing_percent = 100 * p)
}, error = function(e) tibble(error = conditionMessage(e)))
write_csv(pred_primary, file.path(OUT_DIR, "missingness_model_primary_predicted_activity_probabilities.csv"))

# ----------------------------
# Predicted probabilities for activity3 from all-epochs sensitivity model
# ----------------------------
new_all <- model_all_dat %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
  slice(rep(1, length(levels(model_all_dat$activity3))))
new_all$activity3 <- factor(levels(model_all_dat$activity3), levels = levels(model_all_dat$activity3))
new_all$p_id <- model_all_dat$p_id[1]
for (v in c("day_period", "day_type", "weather_info", "gender")) {
  if (v %in% names(model_all_dat)) {
    ref <- levels(model_all_dat[[v]])[1]
    new_all[[v]] <- factor(ref, levels = levels(model_all_dat[[v]]))
  }
}

pred_all <- tryCatch({
  if (!fit_all$fallback && HAS_LME4) {
    p <- predict(fit_all$model, newdata = new_all, type = "response", re.form = NA, allow.new.levels = TRUE)
  } else {
    p <- predict(fit_all$model, newdata = new_all, type = "response")
  }
  tibble(activity3 = new_all$activity3, predicted_missing_rate = p, predicted_missing_percent = 100 * p)
}, error = function(e) tibble(error = conditionMessage(e)))
write_csv(pred_all, file.path(OUT_DIR, "missingness_model_all_epochs_predicted_activity_probabilities.csv"))

# ----------------------------
# Figures
# ----------------------------
plot_activity <- activity_summary %>%
  mutate(activity3 = forcats::fct_reorder(activity3, missing_percent)) %>%
  ggplot(aes(x = activity3, y = missing_percent)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  coord_flip() +
  labs(
    title = "Raw heart-rate missingness by activity stratum",
    subtitle = "Error bars show approximate 95% confidence intervals for row-level missingness rates.",
    x = NULL,
    y = "Missing raw HR (%)"
  ) +
  theme_minimal(base_size = 12)

plot_primary_forest <- fit_primary$or_table %>%
  filter(!str_detect(term, "^p_id")) %>%
  mutate(term_label = clean_term_labels(term), term_label = forcats::fct_reorder(term_label, odds_ratio)) %>%
  ggplot(aes(x = odds_ratio, y = term_label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Supplementary primary model: driving vs non-driving sedentary",
    subtitle = "Odds ratios for missing raw HR. Reference is non-driving sedentary.",
    x = "Odds ratio for missing raw HR (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10)

plot_all_forest <- fit_all$or_table %>%
  filter(!str_detect(term, "^p_id")) %>%
  mutate(term_label = clean_term_labels(term), term_label = forcats::fct_reorder(term_label, odds_ratio)) %>%
  ggplot(aes(x = odds_ratio, y = term_label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "All-epochs sensitivity missingness model",
    subtitle = "Odds ratios for missing raw HR. Reference activity is non-driving sedentary.",
    x = "Odds ratio for missing raw HR (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10)

plot_driving_intensity_forest <- fit_driving_intensity$or_table %>%
  filter(!str_detect(term, "^p_id")) %>%
  mutate(term_label = clean_term_labels(term), term_label = forcats::fct_reorder(term_label, odds_ratio)) %>%
  ggplot(aes(x = odds_ratio, y = term_label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Supplementary driving-only missingness model",
    subtitle = "Driving epochs only. Numeric predictors are per 1 SD.",
    x = "Odds ratio for missing raw HR (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10)

plot_driving_forest <- fit_driving$or_table %>%
  filter(!str_detect(term, "^p_id")) %>%
  mutate(term_label = clean_term_labels(term), term_label = forcats::fct_reorder(term_label, odds_ratio)) %>%
  ggplot(aes(x = odds_ratio, y = term_label)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    title = "Adjusted driving-only sensitivity missingness model",
    subtitle = "Driving epochs only; includes weather, time context, and participant traits. Numeric variables are per 1 SD.",
    x = "Odds ratio for missing raw HR (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10)

plot_context <- context_summary %>%
  filter(context %in% c("weather_info", "speed_quartile", "jam_factor_quartile", "acceleration_energy_quartile", "rotational_energy_quartile")) %>%
  mutate(level = forcats::fct_inorder(level)) %>%
  ggplot(aes(x = level, y = missing_percent)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  facet_wrap(~ context, scales = "free_x", ncol = 1) +
  labs(
    title = "Observed driving HR missingness across contextual bins",
    subtitle = "Quartiles are computed within driving epochs.",
    x = NULL,
    y = "Missing raw HR (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

pdf_file <- file.path(FIG_DIR, "Figure_Supplementary_HR_Missingness_Diagnostics.pdf")
png_file <- file.path(FIG_DIR, "Figure_Supplementary_HR_Missingness_Diagnostics.png")

grDevices::pdf(pdf_file, width = 8.5, height = 6.5)
print(plot_activity)
print(plot_primary_forest)
print(plot_all_forest)
print(plot_driving_intensity_forest)
print(plot_driving_forest)
print(plot_context)
grDevices::dev.off()

grDevices::png(png_file, width = 2200, height = 1600, res = 220)
print(plot_activity)
print(plot_primary_forest)
print(plot_driving_intensity_forest)
print(plot_context)
grDevices::dev.off()

# ----------------------------
# Human-readable notes
# ----------------------------
notes <- c(
  "HR missingness analysis",
  "=======================",
  paste0("Input file: ", DATA_FILE),
  paste0("Output folder: ", OUT_DIR),
  "",
  "Outcome:",
  "  hr_missing = 1 if raw_hr is missing; 0 otherwise.",
  "",
  "Supplementary descriptive comparison:",
  "  missingness_driving_vs_sedentary_summary.csv reports row-level missingness for the driving versus non-driving sedentary contrast reported in the Supplement.",
  "  missingness_activity_summary.csv reports row-level missingness by all activity3 strata.",
  "",
  "Supplementary primary model:",
  paste0("  ", deparse(form_primary_mixed)),
  "  Reference activity level is non_driving_sedentary.",
  "",
  "All-epochs sensitivity model:",
  paste0("  ", deparse(form_all_mixed)),
  "  Reference activity level is non_driving_sedentary.",
  "",
  "Row-retention audit:",
  "  missingness_model_row_retention_audit.csv reports rows, subjects, and missingness rates before and after complete-case filtering for each model.",
  "",
  "Supplementary driving-only model:",
  paste0("  ", deparse(form_driving_intensity_mixed)),
  "  Numeric driving-intensity predictors are z-scored; ORs are per 1 SD increment.",
  "  This model corresponds to the driving-only analysis summarized in Supplementary Table S1.",
  "",
  "Adjusted driving-only context sensitivity model:",
  paste0("  ", deparse(form_driving_mixed)),
  "  Numeric driving context predictors are z-scored; ORs are per 1 SD increment.",
  "",
  "Reproducibility notes:",
  "  The primary and driving-only models correspond to the two inferential analyses summarized in Supplementary Table S1.",
  "  The all-activity and adjusted driving-only models are retained as sensitivity analyses and should not replace the reported Table S1 estimates.",
  "  Model-specific complete-case retention is documented in missingness_model_row_retention_audit.csv."
)
writeLines(notes, file.path(OUT_DIR, "missingness_analysis_notes.txt"))

message("Done. Results written to: ", OUT_DIR)
message("Main figure: ", pdf_file)
