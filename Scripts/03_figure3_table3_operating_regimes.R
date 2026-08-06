# ============================================================
# 03_figure3_table3_operating_regimes.R
#
# PURPOSE
#   Reproduce Figure 3 and Table 3 for the revised
#   npj Digital Public Health manuscript:
#
#     Wearable sensing reveals the structure of cardiac
#     activation associated with everyday driving
#
#   Figure 3 compares baseline-referenced cardiac operating regimes
#   during DRIVING and NONDRIVING_SEDENTARY, overall and separately
#   for weekdays and weekends.
#
#   Panels:
#     A. Distribution of 60-s normalized heart rate (NHR)
#     B. Participant-level percentage of observed time with NHR > 0
#     C. Participant-level standard deviation of NHR
#     D. Participant-level returns to baseline per observed hour
#
#   Table 3 reports participant-level means with 5th-95th percentile
#   ranges and paired t-test p values for the same three calendar strata.
#
# ANALYTIC DEFINITIONS
#   Participant-level baseline:
#
#     HR_base,i = mean_d(HR_base,id)
#
#   Normalized heart rate:
#
#     NHR_itd = HR_raw,itd - HR_base,i
#
#   A return to baseline is a transition from NHR > 0 to NHR <= 0
#   between consecutive available observations within the same participant,
#   behavioral context, and calendar stratum. Rates use observed hours as
#   the denominator.
#
# INPUT
#   Data/NUBI_Data_<resolution>sec_Level_MASTER_CLEAN.csv
#
# DEFAULT
#   NUBI_RES_SECONDS=60
#
# MAJOR OUTPUTS
#   Results/paper_figs/<timestamp>_<resolution>sec_figure3_table3_operating_regimes/
#
#     Figures/Figure3_Operating_Regimes.pdf
#     Figures/Figure3_Operating_Regimes.png
#     Tables/Table3_Operating_Regimes.csv
#     Tables/Table3_Operating_Regimes.tex
#     Diagnostics/...
#     run_log.txt
#
# REPOSITORY SCOPE
#   The public repository starts from curated analysis-ready datasets and
#   does not reconstruct them from raw wearable, smartphone, vehicle, or
#   ground-truth streams.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

options(warn = 1)
# ----------------------------
# User-configurable settings
# ----------------------------
LOCAL_TZ <- "America/Chicago"

RES_SECONDS <- suppressWarnings(
  as.integer(Sys.getenv("NUBI_RES_SECONDS", unset = "60"))
)
if (is.na(RES_SECONDS) || !RES_SECONDS %in% c(10L, 30L, 60L)) {
  stop("NUBI_RES_SECONDS must be one of 10, 30, or 60.")
}

ACT3_DRIVING <- "driving"
ACT3_ND_SED  <- "non_driving_sedentary"

USE_COMMON_SUBJECTS_FOR_BOTH_ACTIVITIES <- TRUE
MIN_ROWS_PER_SUBJECT_CELL <- 2L
DENSITY_TRIM_Q <- c(0.001, 0.999)

RETURN_GAP_TOLERANCE_SECONDS <- as.numeric(RES_SECONDS) * 1.5

COL_DRIVING <- "#E69F00"
COL_SED     <- "#7F7F7F"
COL_ERRBAR  <- "#2F2F2F"
BAR_ALPHA <- 0.82

PAL_ACTIVITY <- c(
  "driving" = COL_DRIVING,
  "sedentary" = COL_SED
)

ACTIVITY_LABELS <- c(
  "sedentary" = "Non-driving
sedentary",
  "driving" = "Driving"
)

STRATUM_LEVELS <- c("Overall", "Weekdays", "Weekends")

# ----------------------------
# Resolve repository paths
# ----------------------------
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)

this_script <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
} else {
  tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = TRUE),
    error = function(e) NA_character_
  )
}

if (is.na(this_script)) {
  stop(
    "Could not determine the script location. ",
    "Run with Rscript or source the file from the repository Scripts/ folder."
  )
}

script_dir <- dirname(this_script)
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

message("Script directory: ", script_dir)
message("Repository root: ", project_root)

in_path <- switch(
  as.character(RES_SECONDS),
  "10" = file.path(project_root, "Data", "NUBI_Data_10sec_Level_MASTER_CLEAN.csv"),
  "30" = file.path(project_root, "Data", "NUBI_Data_30sec_Level_MASTER_CLEAN.csv"),
  "60" = file.path(project_root, "Data", "NUBI_Data_60sec_Level_MASTER_CLEAN.csv")
)

if (!file.exists(in_path)) {
  stop(
    "Dataset not found: ", in_path, "\n",
    "The public repository currently includes the 60-second curated dataset. ",
    "Use 60 seconds, or add the corresponding ", RES_SECONDS,
    "-second dataset under Data/."
  )
}

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(
  project_root,
  "Results", "paper_figs",
  paste0(stamp, "_", RES_SECONDS, "sec_figure3_table3_operating_regimes")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "Figures")
table_dir <- file.path(out_dir, "Tables")
diag_dir <- file.path(out_dir, "Diagnostics")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "run_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_save_pdf <- function(plot_obj, path, w = 7, h = 5) {
  ok <- tryCatch({
    ggsave(filename = path, plot = plot_obj, width = w, height = h, device = grDevices::cairo_pdf)
    TRUE
  }, error = function(e) FALSE)
  if (!ok) ggsave(filename = path, plot = plot_obj, width = w, height = h, device = "pdf", useDingbats = FALSE)
}

safe_save_png <- function(plot_obj, path, w = 7, h = 5, dpi = 300) {
  ggsave(filename = path, plot = plot_obj, width = w, height = h, dpi = dpi)
}

log_msg("SCRIPT VERSION: GITHUB-READY MANUSCRIPT-SYNCHRONIZED FIXED 2026-07-23")
log_msg("Script: 03_figure3_table3_operating_regimes.R")
log_msg("Starting Figure 3 and Table 3 operating-regime analysis")
log_msg("Input: ", in_path)
log_msg("Output: ", out_dir)
log_msg("USE_COMMON_SUBJECTS_FOR_BOTH_ACTIVITIES: ", USE_COMMON_SUBJECTS_FOR_BOTH_ACTIVITIES)

# ============================================================
# HELPERS
# ============================================================
snakeify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

parse_time_local <- function(x, tz = LOCAL_TZ) {
  # Match Figure 2 convention: strip trailing Z and interpret displayed clock
  # time as local America/Chicago wall-clock time.
  if (inherits(x, "POSIXt")) return(lubridate::force_tz(x, tzone = tz))
  x <- as.character(x)
  x2 <- trimws(sub("[Zz]$", "", x))
  tt <- suppressWarnings(as.POSIXct(x2, format = "%Y-%m-%dT%H:%M:%S", tz = tz))
  if (sum(!is.na(tt)) > 0) return(tt)
  tt <- suppressWarnings(as.POSIXct(x2, format = "%Y-%m-%d %H:%M:%S", tz = tz))
  if (sum(!is.na(tt)) > 0) return(tt)
  suppressWarnings(lubridate::ymd_hms(x2, tz = tz, quiet = TRUE))
}

canonicalize_names <- function(dt) {
  stopifnot(is.data.table(dt))
  old <- names(dt)
  sn <- snakeify(old)
  canon <- c(
    p_id = "p_id", pid = "p_id", participant_id = "p_id", participantid = "p_id",
    time = "time", timestamp = "time", datetime = "time", date_time = "time",
    raw_hr = "raw_hr", hr = "raw_hr",
    bl_hr = "bl_hr", hr_bl = "bl_hr", hrbl = "bl_hr",
    baseline_hr = "bl_hr", hr_baseline = "bl_hr",
    activity3 = "activity3", activity = "activity",
    day_num = "day_num", daynum = "day_num", days = "days",
    day_type = "day_type"
  )
  new <- sn
  hit <- sn %in% names(canon)
  new[hit] <- unname(canon[sn[hit]])
  if (!identical(old, new)) setnames(dt, old, make.unique(new, sep = "_"))
  dt
}

numify <- function(x) suppressWarnings(as.numeric(x))

make_day_key <- function(dt, tz = LOCAL_TZ) {
  if ("day_num" %in% names(dt)) {
    x <- dt[["day_num"]]
    if (is.numeric(x) || is.integer(x)) return(as.integer(x))
    dig <- suppressWarnings(as.integer(stringr::str_extract(as.character(x), "\\d+")))
    if (!all(is.na(dig))) return(dig)
  }
  if ("days" %in% names(dt)) {
    x <- dt[["days"]]
    if (is.numeric(x) || is.integer(x)) return(as.integer(x))
    dig <- suppressWarnings(as.integer(stringr::str_extract(as.character(x), "\\d+")))
    if (!all(is.na(dig))) return(dig)
  }
  if ("dt_time" %in% names(dt) && inherits(dt[["dt_time"]], "POSIXt")) {
    d <- as.Date(dt[["dt_time"]], tz = tz)
    u <- sort(unique(d))
    return(as.integer(match(d, u)))
  }
  NULL
}

agg_baseline <- function(x, fun = "median") {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (fun == "median") return(stats::median(x, na.rm = TRUE))
  if (fun == "mean") return(mean(x, na.rm = TRUE))
  stop("Unsupported baseline function: ", fun)
}

fmt_sig <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "fg", flag = "#"))
}

# Fixed-decimal formatter for publication tables. This prevents values such as
# 100. or 0 from appearing in LaTeX output and keeps each metric internally
# consistent across means and percentile limits.
fmt_fixed <- function(x, digits = 1) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "f"))
}

fmt_p <- function(p, digits = 3) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "<0.001", formatC(p, digits = digits, format = "fg", flag = "#"))
  )
}

fmt_ci <- function(lower, upper, digits = 3) {
  ifelse(
    is.na(lower) | is.na(upper),
    NA_character_,
    paste0("[", fmt_sig(lower, digits), ", ", fmt_sig(upper, digits), "]")
  )
}

round_table_3digits <- function(df) {
  out <- as.data.frame(df)
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      if (tolower(nm) %in% c("p.value", "p_value", "p")) {
        out[[nm]] <- fmt_p(out[[nm]], digits = 3)
      } else {
        out[[nm]] <- fmt_sig(out[[nm]], digits = 3)
      }
    }
  }
  out
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x
}

write_simple_latex_table <- function(df, path, caption, label) {
  df <- as.data.frame(df)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("\\begin{table}[!ht]\n\\centering\n", file = con)
  cat("\\caption{", caption, "}\n", sep = "", file = con)
  cat("\\label{", label, "}\n", sep = "", file = con)
  cat("\\begin{tabular}{", paste(rep("l", ncol(df)), collapse = ""), "}\n", sep = "", file = con)
  cat("\\toprule\n", file = con)
  cat(paste(latex_escape(names(df)), collapse = " & "), " \\\\\n", sep = "", file = con)
  cat("\\midrule\n", file = con)
  for (i in seq_len(nrow(df))) {
    cat(paste(latex_escape(df[i, , drop = TRUE]), collapse = " & "), " \\\\\n", sep = "", file = con)
  }
  cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n", file = con)
}

mean_ci <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  se <- stats::sd(x) / sqrt(n)
  tcrit <- stats::qt(0.975, df = n - 1)
  c(n = n, mean = m, se = se, lower = m - tcrit * se, upper = m + tcrit * se)
}

# ============================================================
# READ + STANDARDIZE
# ============================================================
log_msg("Reading data ...")
dt <- fread(in_path, showProgress = TRUE)
dt <- canonicalize_names(dt)

need <- c("p_id", "time", "raw_hr", "bl_hr", "activity3")
miss <- setdiff(need, names(dt))
if (length(miss) > 0) stop("Missing required columns after standardization: ", paste(miss, collapse = ", "))

dt[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]
dt <- dt[!is.na(dt_time)]
dt[, raw_hr_num := numify(raw_hr)]
dt[, bl_hr_num := numify(bl_hr)]
dt[, activity3 := trimws(tolower(as.character(activity3)))]
dt[, day_num_num := make_day_key(.SD, tz = LOCAL_TZ)]

if (all(is.na(dt$day_num_num))) stop("Could not construct day key from day_num/days/time.")

# ============================================================
# PERSON-LEVEL BASELINE WITH PARTICIPANT-DAY QC
# ============================================================
# Match the predictive decomposition / ML scripts in substance:
#   1) recover the single Apple Watch/iPhone baseline value for each
#      participant x day;
#   2) compute one participant-level baseline as the mean of those daily
#      baseline values.
#
# Apple Watch/iPhone should provide one baseline HR per participant-day, repeated
# across rows from that day. We therefore check that the repeated value is
# actually constant within each participant x day and then take that unique
# value directly, rather than using a median as a fallback. This adds explicit
# quality control to the baseline construction.
BL_HR_DAY_UNIQUE_TOL <- 1e-8

baseline_by_subj_day_qc <- dt[
  is.finite(bl_hr_num) & !is.na(day_num_num),
  .(
    n_rows_used = .N,
    n_unique_bl_hr = uniqueN(bl_hr_num),
    bl_hr_min = min(bl_hr_num),
    bl_hr_max = max(bl_hr_num),
    bl_hr_range = max(bl_hr_num) - min(bl_hr_num),
    bl_hr_values = paste(sort(unique(bl_hr_num)), collapse = "; ")
  ),
  by = .(p_id, day_num = as.integer(day_num_num))
]

write_csv(
  as.data.frame(baseline_by_subj_day_qc),
  file.path(out_dir, "expanded_figure3_baseline_by_subject_day_QC.csv")
)
write_csv(
  round_table_3digits(baseline_by_subj_day_qc),
  file.path(out_dir, "expanded_figure3_baseline_by_subject_day_QC_3digits.csv")
)

if (nrow(baseline_by_subj_day_qc) == 0) {
  stop("No finite baseline values available to compute participant-level baselines.")
}

bad_bl_days <- baseline_by_subj_day_qc[bl_hr_range > BL_HR_DAY_UNIQUE_TOL]
if (nrow(bad_bl_days) > 0) {
  write_csv(
    as.data.frame(bad_bl_days),
    file.path(out_dir, "expanded_figure3_baseline_by_subject_day_QC_FAILED.csv")
  )
  stop(
    "Baseline QC failed: at least one participant-day has more than one finite bl_hr value. ",
    "See expanded_figure3_baseline_by_subject_day_QC_FAILED.csv in the output directory."
  )
}

baseline_by_subj_day <- baseline_by_subj_day_qc[
  ,
  .(
    bl_hr_day = bl_hr_min,
    n_rows_used = n_rows_used,
    n_unique_bl_hr = n_unique_bl_hr,
    bl_hr_range = bl_hr_range
  ),
  by = .(p_id, day_num)
]

write_csv(as.data.frame(baseline_by_subj_day), file.path(out_dir, "expanded_figure3_baseline_by_subject_day.csv"))
write_csv(round_table_3digits(baseline_by_subj_day), file.path(out_dir, "expanded_figure3_baseline_by_subject_day_3digits.csv"))

subj_bl <- baseline_by_subj_day[
  ,
  .(
    bl_hr_subject = mean(bl_hr_day, na.rm = TRUE),
    n_days_with_bl = .N,
    bl_hr_day_sd = stats::sd(bl_hr_day, na.rm = TRUE)
  ),
  by = .(p_id)
]
write_csv(as.data.frame(subj_bl), file.path(out_dir, "expanded_figure3_subject_baseline_summary.csv"))
write_csv(round_table_3digits(subj_bl), file.path(out_dir, "expanded_figure3_subject_baseline_summary_3digits.csv"))
log_msg("Wrote participant-level baseline summary: n_subjects=", nrow(subj_bl))

dt <- dt[activity3 %in% c(ACT3_DRIVING, ACT3_ND_SED)]
if (nrow(dt) == 0) stop("No rows remain after filtering to driving and non-driving sedentary.")

log_msg("Rows after activity filter: ", nrow(dt))
log_msg("Subjects after activity filter: ", uniqueN(dt$p_id))
capture.output(print(dt[, .N, by = activity3][order(activity3)]), file = log_file, append = TRUE)

dd0 <- dt[
  is.finite(raw_hr_num) & !is.na(day_num_num),
  .(
    p_id,
    time = as.character(time),
    dt_time,
    raw_hr = raw_hr_num,
    bl_hr = bl_hr_num,
    activity3,
    day_num = as.integer(day_num_num),
    day_type = if ("day_type" %in% names(dt)) as.character(day_type) else NA_character_
  )
]
if (nrow(dd0) == 0) stop("No complete rows with finite raw_hr and day_num.")

# Attach the participant-level baseline computed above and calculate NHR.
# This matches the predictive decomposition definition in the paper:
# NHR_itd = HR_raw,itd - HR_base,i.
dd <- merge(dd0, subj_bl[, .(p_id, bl_hr_subject, n_days_with_bl)], by = "p_id", all.x = TRUE)
dd <- dd[is.finite(bl_hr_subject)]
dd[, nhr := raw_hr - bl_hr_subject]
dd[, nhr_gt_zero := nhr > 0]

# Study-day labels follow the manuscript's Monday-Sunday convention:
# days 1-5 are weekdays and days 6-7 are weekends.
if (any(!dd$day_num %in% 1:7)) {
  stop("Study-day values outside 1-7 were detected after preprocessing.")
}
dd[, day_type_analysis := fifelse(day_num %in% 1:5, "weekdays", "weekend")]

dd[, activity_binary := fifelse(activity3 == ACT3_DRIVING, "driving", "sedentary")]
dd[, activity_binary := factor(activity_binary, levels = c("sedentary", "driving"))]
dd[, day_type_analysis := factor(day_type_analysis, levels = c("weekdays", "weekend"))]

if (USE_COMMON_SUBJECTS_FOR_BOTH_ACTIVITIES) {
  subj_drive <- unique(dd[activity3 == ACT3_DRIVING, p_id])
  subj_sed <- unique(dd[activity3 == ACT3_ND_SED, p_id])
  common_subj <- intersect(subj_drive, subj_sed)
  dd <- dd[p_id %in% common_subj]
  log_msg("Common-subject restriction for both activities applied.")
  log_msg("Common subjects: ", length(common_subj))
}

log_msg("Analysis rows: ", nrow(dd))
log_msg("Analysis subjects: ", uniqueN(dd$p_id))
capture.output(print(dd[, .(n_rows = .N, n_subjects = uniqueN(p_id)), by = .(activity_binary, day_type_analysis)]), file = log_file, append = TRUE)

fwrite(dd, file.path(out_dir, "expanded_figure3_analysis_data.csv"))

# ============================================================
# CREATE OVERALL + WEEKDAY + WEEKEND DATASETS
# ============================================================
dd_week <- copy(dd)
dd_week[, stratum := fifelse(day_type_analysis == "weekdays", "Weekdays", "Weekends")]

dd_overall <- copy(dd)
dd_overall[, stratum := "Overall"]

dd_exp <- rbindlist(list(dd_overall, dd_week), use.names = TRUE, fill = TRUE)
dd_exp[, stratum := factor(stratum, levels = STRATUM_LEVELS)]
dd_exp[, activity_binary := factor(activity_binary, levels = c("sedentary", "driving"))]
dd_exp[, activity_label := factor(unname(ACTIVITY_LABELS[as.character(activity_binary)]), levels = unname(ACTIVITY_LABELS))]

# ============================================================
# RETURNS TO BASELINE
# ============================================================
# Count only downward crossings from above baseline to at/below baseline.
# The denominator is observed time, not elapsed clock time, so the rate is
# events per hour of available 60-s observations within each cell.
setorder(dd_exp, p_id, stratum, activity_binary, dt_time)
dd_exp[, `:=`(
  prev_nhr_gt_zero = shift(nhr_gt_zero),
  prev_time = shift(dt_time),
  gap_seconds = as.numeric(difftime(dt_time, shift(dt_time), units = "secs"))
), by = .(p_id, stratum, activity_binary)]

dd_exp[, return_to_baseline :=
  !is.na(prev_nhr_gt_zero) & prev_nhr_gt_zero & !nhr_gt_zero &
  is.finite(gap_seconds) & gap_seconds > 0 &
  gap_seconds <= RETURN_GAP_TOLERANCE_SECONDS]

return_events <- dd_exp[return_to_baseline == TRUE,
  .(p_id, stratum, activity_binary, dt_time, nhr, gap_seconds)]
fwrite(return_events, file.path(out_dir, "expanded_figure3_return_events.csv"))

# ============================================================
# TABLES FOR FIGURE 3
# ============================================================
row_summary <- dd_exp %>%
  as.data.frame() %>%
  group_by(stratum, activity_binary) %>%
  summarise(
    n_rows = n(),
    n_subjects = n_distinct(p_id),
    n_subject_days = n_distinct(paste(p_id, day_num, sep = "__")),
    mean_nhr = mean(nhr, na.rm = TRUE),
    median_nhr = median(nhr, na.rm = TRUE),
    sd_nhr_rows = sd(nhr, na.rm = TRUE),
    pct_nhr_gt_zero_rows = mean(nhr_gt_zero, na.rm = TRUE),
    q05_nhr = as.numeric(quantile(nhr, 0.05, na.rm = TRUE)),
    q25_nhr = as.numeric(quantile(nhr, 0.25, na.rm = TRUE)),
    q50_nhr = as.numeric(quantile(nhr, 0.50, na.rm = TRUE)),
    q75_nhr = as.numeric(quantile(nhr, 0.75, na.rm = TRUE)),
    q95_nhr = as.numeric(quantile(nhr, 0.95, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(stratum, activity_binary)

write_csv(row_summary, file.path(out_dir, "expanded_figure3_row_level_summary.csv"))
write_csv(round_table_3digits(row_summary), file.path(out_dir, "expanded_figure3_row_level_summary_3digits.csv"))

subject_cell <- dd_exp[
  ,
  .(
    n_rows = .N,
    n_above_zero = sum(nhr_gt_zero, na.rm = TRUE),
    pct_nhr_gt_zero = mean(nhr_gt_zero, na.rm = TRUE),
    mean_nhr = mean(nhr, na.rm = TRUE),
    median_nhr = median(nhr, na.rm = TRUE),
    sd_nhr = sd(nhr, na.rm = TRUE),
    n_return_to_baseline = sum(return_to_baseline, na.rm = TRUE),
    observed_hours = .N * RES_SECONDS / 3600,
    # data.table does not expose newly created summary aliases inside the same
    # .() expression, so compute the rate directly from source quantities.
    returns_per_hour = sum(return_to_baseline, na.rm = TRUE) / (.N * RES_SECONDS / 3600),
    q05_nhr = as.numeric(quantile(nhr, 0.05, na.rm = TRUE)),
    q50_nhr = as.numeric(quantile(nhr, 0.50, na.rm = TRUE)),
    q95_nhr = as.numeric(quantile(nhr, 0.95, na.rm = TRUE))
  ),
  by = .(p_id, stratum, activity_binary)
]
subject_cell <- subject_cell[n_rows >= MIN_ROWS_PER_SUBJECT_CELL & is.finite(sd_nhr)]
subject_cell[, activity_label := factor(unname(ACTIVITY_LABELS[as.character(activity_binary)]), levels = unname(ACTIVITY_LABELS))]
subject_cell[, stratum := factor(stratum, levels = STRATUM_LEVELS)]

write_csv(as.data.frame(subject_cell), file.path(out_dir, "expanded_figure3_subject_cell_summary.csv"))
write_csv(round_table_3digits(subject_cell), file.path(out_dir, "expanded_figure3_subject_cell_summary_3digits.csv"))

panel_values <- subject_cell %>%
  as.data.frame() %>%
  group_by(stratum, activity_binary) %>%
  summarise(
    n_subjects = n_distinct(p_id),
    mean_subject_pct_nhr_gt_zero = mean(pct_nhr_gt_zero, na.rm = TRUE),
    q05_subject_pct_nhr_gt_zero = as.numeric(quantile(pct_nhr_gt_zero, 0.05, na.rm = TRUE)),
    q50_subject_pct_nhr_gt_zero = as.numeric(quantile(pct_nhr_gt_zero, 0.50, na.rm = TRUE)),
    q95_subject_pct_nhr_gt_zero = as.numeric(quantile(pct_nhr_gt_zero, 0.95, na.rm = TRUE)),
    mean_subject_sd_nhr = mean(sd_nhr, na.rm = TRUE),
    q05_subject_sd_nhr = as.numeric(quantile(sd_nhr, 0.05, na.rm = TRUE)),
    q50_subject_sd_nhr = as.numeric(quantile(sd_nhr, 0.50, na.rm = TRUE)),
    q95_subject_sd_nhr = as.numeric(quantile(sd_nhr, 0.95, na.rm = TRUE)),
    mean_subject_mean_nhr = mean(mean_nhr, na.rm = TRUE),
    q05_subject_mean_nhr = as.numeric(quantile(mean_nhr, 0.05, na.rm = TRUE)),
    q50_subject_mean_nhr = as.numeric(quantile(mean_nhr, 0.50, na.rm = TRUE)),
    q95_subject_mean_nhr = as.numeric(quantile(mean_nhr, 0.95, na.rm = TRUE)),
    mean_subject_returns_per_hour = mean(returns_per_hour, na.rm = TRUE),
    q05_subject_returns_per_hour = as.numeric(quantile(returns_per_hour, 0.05, na.rm = TRUE)),
    q50_subject_returns_per_hour = as.numeric(quantile(returns_per_hour, 0.50, na.rm = TRUE)),
    q95_subject_returns_per_hour = as.numeric(quantile(returns_per_hour, 0.95, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(stratum, activity_binary)

write_csv(panel_values, file.path(out_dir, "expanded_figure3_panel_values.csv"))
write_csv(round_table_3digits(panel_values), file.path(out_dir, "expanded_figure3_panel_values_3digits.csv"))

panel_values_paper <- panel_values %>%
  mutate(
    Activity = recode(as.character(activity_binary), sedentary = "Non-driving sedentary", driving = "Driving"),
    Stratum = as.character(stratum),
    `% time NHR > 0` = paste0(fmt_sig(100 * mean_subject_pct_nhr_gt_zero, 3), "%"),
    `% time NHR > 0, 5th-95th` = paste0("[", fmt_sig(100 * q05_subject_pct_nhr_gt_zero, 3), "%, ", fmt_sig(100 * q95_subject_pct_nhr_gt_zero, 3), "%]"),
    `SD(NHR) [bpm]` = fmt_sig(mean_subject_sd_nhr, 3),
    `SD(NHR), 5th-95th [bpm]` = fmt_ci(q05_subject_sd_nhr, q95_subject_sd_nhr, 3),
    `Mean NHR [bpm]` = fmt_sig(mean_subject_mean_nhr, 3),
    `Mean NHR, 5th-95th [bpm]` = fmt_ci(q05_subject_mean_nhr, q95_subject_mean_nhr, 3),
    `Returns to baseline [events/h]` = fmt_sig(mean_subject_returns_per_hour, 3),
    `Returns to baseline, 5th-95th [events/h]` = fmt_ci(q05_subject_returns_per_hour, q95_subject_returns_per_hour, 3)
  ) %>%
  select(
    Stratum,
    Activity,
    n_subjects,
    `% time NHR > 0`,
    `% time NHR > 0, 5th-95th`,
    `SD(NHR) [bpm]`,
    `SD(NHR), 5th-95th [bpm]`,
    `Mean NHR [bpm]`,
    `Mean NHR, 5th-95th [bpm]`,
    `Returns to baseline [events/h]`,
    `Returns to baseline, 5th-95th [events/h]`
  )

write_csv(panel_values_paper, file.path(out_dir, "expanded_figure3_panel_values_paper_table_3digits.csv"))
write_simple_latex_table(
  panel_values_paper,
  file.path(out_dir, "expanded_figure3_latex_panel_values_3digits.tex"),
  caption = "Expanded Figure 3 participant-level summaries by activity and calendar stratum.",
  label = "tab:expanded_figure3_panel_values"
)

# Paired participant-level driving - sedentary contrasts within each stratum.
paired_contrasts <- subject_cell %>%
  as.data.frame() %>%
  select(p_id, stratum, activity_binary, pct_nhr_gt_zero, sd_nhr, mean_nhr, returns_per_hour) %>%
  pivot_wider(
    names_from = activity_binary,
    values_from = c(pct_nhr_gt_zero, sd_nhr, mean_nhr, returns_per_hour)
  ) %>%
  filter(
    is.finite(pct_nhr_gt_zero_driving), is.finite(pct_nhr_gt_zero_sedentary),
    is.finite(sd_nhr_driving), is.finite(sd_nhr_sedentary),
    is.finite(mean_nhr_driving), is.finite(mean_nhr_sedentary),
    is.finite(returns_per_hour_driving), is.finite(returns_per_hour_sedentary)
  ) %>%
  mutate(
    diff_pct_nhr_gt_zero = pct_nhr_gt_zero_driving - pct_nhr_gt_zero_sedentary,
    diff_sd_nhr = sd_nhr_driving - sd_nhr_sedentary,
    diff_mean_nhr = mean_nhr_driving - mean_nhr_sedentary,
    diff_returns_per_hour = returns_per_hour_driving - returns_per_hour_sedentary
  )

contrast_stats <- paired_contrasts %>%
  group_by(stratum) %>%
  summarise(
    n_paired_subjects = n(),
    pct_diff_mean = mean(diff_pct_nhr_gt_zero, na.rm = TRUE),
    pct_diff_se = sd(diff_pct_nhr_gt_zero, na.rm = TRUE) / sqrt(n()),
    pct_diff_lower = pct_diff_mean - qt(0.975, df = n() - 1) * pct_diff_se,
    pct_diff_upper = pct_diff_mean + qt(0.975, df = n() - 1) * pct_diff_se,
    pct_diff_p = t.test(diff_pct_nhr_gt_zero)$p.value,
    sd_diff_mean = mean(diff_sd_nhr, na.rm = TRUE),
    sd_diff_se = sd(diff_sd_nhr, na.rm = TRUE) / sqrt(n()),
    sd_diff_lower = sd_diff_mean - qt(0.975, df = n() - 1) * sd_diff_se,
    sd_diff_upper = sd_diff_mean + qt(0.975, df = n() - 1) * sd_diff_se,
    sd_diff_p = t.test(diff_sd_nhr)$p.value,
    mean_nhr_diff_mean = mean(diff_mean_nhr, na.rm = TRUE),
    mean_nhr_diff_se = sd(diff_mean_nhr, na.rm = TRUE) / sqrt(n()),
    mean_nhr_diff_lower = mean_nhr_diff_mean - qt(0.975, df = n() - 1) * mean_nhr_diff_se,
    mean_nhr_diff_upper = mean_nhr_diff_mean + qt(0.975, df = n() - 1) * mean_nhr_diff_se,
    mean_nhr_diff_p = t.test(diff_mean_nhr)$p.value,
    returns_diff_mean = mean(diff_returns_per_hour, na.rm = TRUE),
    returns_diff_se = sd(diff_returns_per_hour, na.rm = TRUE) / sqrt(n()),
    returns_diff_lower = returns_diff_mean - qt(0.975, df = n() - 1) * returns_diff_se,
    returns_diff_upper = returns_diff_mean + qt(0.975, df = n() - 1) * returns_diff_se,
    returns_diff_p = t.test(diff_returns_per_hour)$p.value,
    .groups = "drop"
  ) %>%
  arrange(stratum)

write_csv(contrast_stats, file.path(out_dir, "expanded_figure3_paired_activity_contrasts.csv"))
write_csv(round_table_3digits(contrast_stats), file.path(out_dir, "expanded_figure3_paired_activity_contrasts_3digits.csv"))

contrast_paper <- contrast_stats %>%
  mutate(
    Stratum = as.character(stratum),
    `n paired subjects` = as.character(n_paired_subjects),
    `Driving - sedentary % time NHR > 0` = paste0(fmt_sig(100 * pct_diff_mean, 3), " pp"),
    `95% CI, percentage points` = paste0("[", fmt_sig(100 * pct_diff_lower, 3), ", ", fmt_sig(100 * pct_diff_upper, 3), "]"),
    `p, % time NHR > 0` = fmt_p(pct_diff_p, 3),
    `Driving - sedentary SD(NHR) [bpm]` = fmt_sig(sd_diff_mean, 3),
    `95% CI, SD(NHR) [bpm]` = fmt_ci(sd_diff_lower, sd_diff_upper, 3),
    `p, SD(NHR)` = fmt_p(sd_diff_p, 3),
    `Driving - sedentary mean NHR [bpm]` = fmt_sig(mean_nhr_diff_mean, 3),
    `95% CI, mean NHR [bpm]` = fmt_ci(mean_nhr_diff_lower, mean_nhr_diff_upper, 3),
    `p, mean NHR` = fmt_p(mean_nhr_diff_p, 3),
    `Driving - sedentary returns [events/h]` = fmt_sig(returns_diff_mean, 3),
    `95% CI, returns [events/h]` = fmt_ci(returns_diff_lower, returns_diff_upper, 3),
    `p, returns` = fmt_p(returns_diff_p, 3)
  ) %>%
  select(
    Stratum,
    `n paired subjects`,
    `Driving - sedentary % time NHR > 0`,
    `95% CI, percentage points`,
    `p, % time NHR > 0`,
    `Driving - sedentary SD(NHR) [bpm]`,
    `95% CI, SD(NHR) [bpm]`,
    `p, SD(NHR)`,
    `Driving - sedentary mean NHR [bpm]`,
    `95% CI, mean NHR [bpm]`,
    `p, mean NHR`,
    `Driving - sedentary returns [events/h]`,
    `95% CI, returns [events/h]`,
    `p, returns`
  )

write_csv(contrast_paper, file.path(out_dir, "expanded_figure3_paired_activity_contrasts_paper_table_3digits.csv"))
write_simple_latex_table(
  contrast_paper,
  file.path(out_dir, "expanded_figure3_latex_paired_contrasts_3digits.tex"),
  caption = "Participant-paired driving-minus-sedentary contrasts for expanded Figure 3.",
  label = "tab:expanded_figure3_paired_contrasts"
)

log_msg("Expanded Figure 3 panel values:")
capture.output(print(panel_values_paper, row.names = FALSE), file = log_file, append = TRUE)
log_msg("Expanded Figure 3 paired contrasts:")
capture.output(print(contrast_paper, row.names = FALSE), file = log_file, append = TRUE)

# ============================================================
# NEW TABLE 3: OVERALL, WEEKDAY, WEEKEND
# ============================================================
# Table 3 now mirrors Figure 3: participant-level means are reported together
# with the 5th--95th percentile range across participants. Because the table
# centers on means, the primary inferential comparison is a paired t-test on
# participant-level values. A paired Wilcoxon signed-rank test is also computed
# and retained in the machine-readable output as a nonparametric robustness
# check. Recovery time is omitted because it is exactly the complement of time
# above baseline.

table3_long <- subject_cell %>%
  as.data.frame() %>%
  select(p_id, stratum, activity_binary, mean_nhr, pct_nhr_gt_zero, sd_nhr, returns_per_hour) %>%
  pivot_longer(
    cols = c(mean_nhr, pct_nhr_gt_zero, sd_nhr, returns_per_hour),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    value = ifelse(metric == "pct_nhr_gt_zero", 100 * value, value),
    metric = factor(
      metric,
      levels = c("mean_nhr", "pct_nhr_gt_zero", "sd_nhr", "returns_per_hour"),
      labels = c(
        "Mean NHR [bpm]",
        "Time with NHR > 0 [%]",
        "SD(NHR) [bpm]",
        "Returns to baseline [events/h]"
      )
    )
  )

table3_summary <- table3_long %>%
  group_by(stratum, metric, activity_binary) %>%
  summarise(
    n_subjects = n_distinct(p_id),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    se = sd / sqrt(n_subjects),
    ci95_lower = mean - qt(0.975, df = n_subjects - 1) * se,
    ci95_upper = mean + qt(0.975, df = n_subjects - 1) * se,
    q05 = as.numeric(quantile(value, 0.05, na.rm = TRUE)),
    q95 = as.numeric(quantile(value, 0.95, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    # The revised manuscript rounds Mean NHR and Time above baseline to one
    # decimal place but displays two decimal positions (for example, 11.20 and
    # 93.70). SD(NHR) and returns to baseline are rounded and displayed to two
    # decimal places. Separate rounding from display precision so the generated
    # table reproduces the manuscript exactly.
    rounding_digits = case_when(
      as.character(metric) %in% c(
        "Mean NHR [bpm]",
        "Time with NHR > 0 [%]"
      ) ~ 1L,
      TRUE ~ 2L
    ),
    display_digits = 2L,
    summary = mapply(
      function(m, lo, hi, round_d, display_d) {
        paste0(
          fmt_fixed(round(m, round_d), display_d),
          " [",
          fmt_fixed(round(lo, round_d), display_d),
          "--",
          fmt_fixed(round(hi, round_d), display_d),
          "]"
        )
      },
      mean, q05, q95, rounding_digits, display_digits,
      USE.NAMES = FALSE
    )
  ) %>%
  select(
    stratum, metric, activity_binary, n_subjects,
    mean, sd, se, ci95_lower, ci95_upper, q05, q95, summary
  )

table3_wide_summary <- table3_summary %>%
  select(stratum, metric, activity_binary, n_subjects, mean, q05, q95, summary) %>%
  pivot_wider(
    names_from = activity_binary,
    values_from = c(n_subjects, mean, q05, q95, summary)
  )

# Paired participant-level tests and effect estimates.
table3_tests <- table3_long %>%
  select(p_id, stratum, metric, activity_binary, value) %>%
  pivot_wider(names_from = activity_binary, values_from = value) %>%
  filter(is.finite(driving), is.finite(sedentary)) %>%
  group_by(stratum, metric) %>%
  summarise(
    n_paired = n(),
    mean_difference = mean(driving - sedentary),
    sd_difference = sd(driving - sedentary),
    se_difference = sd_difference / sqrt(n_paired),
    diff_ci95_lower = mean_difference - qt(0.975, df = n_paired - 1) * se_difference,
    diff_ci95_upper = mean_difference + qt(0.975, df = n_paired - 1) * se_difference,
    paired_t_p = suppressWarnings(t.test(driving, sedentary, paired = TRUE)$p.value),
    paired_wilcoxon_p = suppressWarnings(
      wilcox.test(driving, sedentary, paired = TRUE, exact = FALSE)$p.value
    ),
    .groups = "drop"
  )

table3 <- table3_wide_summary %>%
  left_join(table3_tests, by = c("stratum", "metric")) %>%
  transmute(
    Stratum = as.character(stratum),
    Metric = as.character(metric),
    `Driving, mean [5th--95th percentile]` = summary_driving,
    `Non-driving sedentary, mean [5th--95th percentile]` = summary_sedentary,
    `n paired` = n_paired,
    `Mean paired difference (driving - sedentary)` = mean_difference,
    `95% CI of paired difference, lower` = diff_ci95_lower,
    `95% CI of paired difference, upper` = diff_ci95_upper,
    `Paired t-test p` = paired_t_p,
    `Paired Wilcoxon p` = paired_wilcoxon_p
  ) %>%
  arrange(
    factor(Stratum, levels = STRATUM_LEVELS),
    factor(
      Metric,
      levels = c(
        "Mean NHR [bpm]",
        "Time with NHR > 0 [%]",
        "SD(NHR) [bpm]",
        "Returns to baseline [events/h]"
      )
    )
  )

write_csv(
  table3,
  file.path(table_dir, "Table3_Operating_Regimes.csv")
)
write_csv(
  round_table_3digits(table3),
  file.path(table_dir, "Table3_Operating_Regimes_3digits.csv")
)

# Compact paper-facing table: same means and 5th--95th percentile ranges shown
# in Figure 3, with paired t-test p values. Wilcoxon robustness results and
# paired-difference confidence intervals are retained in machine-readable CSV
# outputs for internal verification and reproducibility.
table3_paper <- table3 %>%
  transmute(
    Stratum,
    Metric,
    Driving = `Driving, mean [5th--95th percentile]`,
    `Non-driving sedentary` = `Non-driving sedentary, mean [5th--95th percentile]`,
    `n` = as.integer(`n paired`),
    `p` = `Paired t-test p`
  )

fmt_p_latex <- function(p, digits = 3) {
  ifelse(
    is.na(p),
    "--",
    ifelse(
      p < 0.001,
      "$<0.001$",
      paste0("$", formatC(p, digits = digits, format = "f"), "$")
    )
  )
}

write_table3_latex <- function(df, path) {
  df <- as.data.frame(df)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)

  cat("\\begin{table*}[!t]\n", file = con)
  cat("\\centering\n", file = con)
  cat(
    "\\caption{Participant-level operating-regime characteristics of baseline-referenced cardiac activation during driving and non-driving sedentary behavior, shown overall and separately for weekdays and weekends. Metrics were calculated separately for each participant from the underlying 60-s observations within each activity and calendar stratum. Values are means [5th--95th percentile range] across participants, matching Figure~\\ref{fig:figure3}; $n$ denotes participants with data in both behavioral contexts, and $p$ values are from paired $t$-tests.}\n",
    file = con
  )
  cat("\\label{tab:operating_regime_overall_weekday_weekend}\n", file = con)
  cat("\\small\n", file = con)
  cat("\\setlength{\\tabcolsep}{4.5pt}\n", file = con)
  cat("\\begin{tabular}{@{}llp{3.25cm}p{3.85cm}cc@{}}\n", file = con)
  cat("\\toprule\n", file = con)
  cat(
    "Stratum & Metric & Driving, mean [5th--95th] & Non-driving sedentary, mean [5th--95th] & $n$ & $p$ \\\\\n",
    file = con
  )
  cat("\\midrule\n", file = con)

  previous_stratum <- NULL
  for (i in seq_len(nrow(df))) {
    current_stratum <- as.character(df$Stratum[i])
    if (!is.null(previous_stratum) && current_stratum != previous_stratum) {
      cat("\\addlinespace\n", file = con)
    }

    metric_tex <- as.character(df$Metric[i])
    metric_tex <- gsub("NHR > 0", "NHR $>$ 0", metric_tex, fixed = TRUE)
    metric_tex <- gsub("[%]", "[\\%]", metric_tex, fixed = TRUE)

    row <- paste0(
      current_stratum, " & ",
      metric_tex, " & ",
      as.character(df$Driving[i]), " & ",
      as.character(df[["Non-driving sedentary"]][i]), " & ",
      as.character(df[["n"]][i]), " & ",
      fmt_p_latex(df[["p"]][i], 3),
      " \\\\"
    )
    cat(row, "\n", file = con, sep = "")
    previous_stratum <- current_stratum
  }

  cat("\\bottomrule\n", file = con)
  cat("\\end{tabular}\n", file = con)
  cat("\\end{table*}\n", file = con)
}

# Verify the published 60-s Table 3 values before writing final outputs.
if (RES_SECONDS == 60L) {
  expected_table3 <- data.frame(
    Stratum = rep(c("Overall", "Weekdays", "Weekends"), each = 4),
    Metric = rep(
      c(
        "Mean NHR [bpm]",
        "Time with NHR > 0 [%]",
        "SD(NHR) [bpm]",
        "Returns to baseline [events/h]"
      ),
      times = 3
    ),
    Driving = c(
      "11.20 [5.90--16.80]",
      "93.70 [81.20--99.60]",
      "7.98 [4.98--12.10]",
      "1.31 [0.22--2.94]",
      "11.10 [5.20--17.60]",
      "93.80 [76.70--100.00]",
      "7.88 [4.41--12.05]",
      "1.34 [0.00--2.96]",
      "12.50 [2.50--22.80]",
      "94.70 [72.90--100.00]",
      "5.99 [2.23--11.91]",
      "1.01 [0.00--4.00]"
    ),
    NonDriving = c(
      "11.50 [4.50--18.80]",
      "81.30 [59.20--94.70]",
      "12.40 [8.98--15.25]",
      "2.46 [1.12--3.96]",
      "11.60 [4.00--19.20]",
      "82.10 [58.90--95.90]",
      "12.11 [8.53--14.59]",
      "2.51 [1.02--4.33]",
      "11.30 [0.70--20.60]",
      "78.40 [40.10--99.20]",
      "11.93 [7.90--17.16]",
      "2.28 [0.21--5.06]"
    ),
    n = c(rep(57L, 8), rep(55L, 4)),
    p_display = c(
      "0.613", "<0.001", "<0.001", "<0.001",
      "0.467", "<0.001", "<0.001", "<0.001",
      "0.313", "<0.001", "<0.001", "<0.001"
    ),
    stringsAsFactors = FALSE
  )

  actual_check <- table3_paper %>%
    transmute(
      Stratum,
      Metric,
      Driving,
      NonDriving = `Non-driving sedentary`,
      n,
      p_display = fmt_p(p, 3)
    ) %>%
    as.data.frame()

  expected_check <- expected_table3
  rownames(actual_check) <- NULL
  rownames(expected_check) <- NULL

  if (!identical(actual_check, expected_check)) {
    write_csv(actual_check, file.path(diag_dir, "Table3_Manuscript_Check_Actual.csv"))
    write_csv(expected_check, file.path(diag_dir, "Table3_Manuscript_Check_Expected.csv"))
    stop(
      "Table 3 manuscript verification failed. ",
      "See Diagnostics/Table3_Manuscript_Check_Actual.csv and ",
      "Diagnostics/Table3_Manuscript_Check_Expected.csv."
    )
  }

  log_msg("Table 3 manuscript checks passed.")
}

write_csv(
  table3_paper %>% mutate(`p` = fmt_p(`p`, 3)),
  file.path(table_dir, "Table3_Operating_Regimes_Manuscript.csv")
)

write_table3_latex(
  table3_paper,
  file.path(table_dir, "Table3_Operating_Regimes.tex")
)

log_msg("New Table 3 (means with 5th--95th percentile ranges):")
capture.output(print(table3_paper, row.names = FALSE), file = log_file, append = TRUE)
log_msg("Table 3 paired-test details, including Wilcoxon robustness checks:")
capture.output(print(table3_tests, row.names = FALSE), file = log_file, append = TRUE)

# ============================================================
# EXPANDED FIGURE 3
# ============================================================
xlims <- quantile(dd_exp$nhr, probs = DENSITY_TRIM_Q, na.rm = TRUE)

density_df <- dd_exp[nhr >= xlims[1] & nhr <= xlims[2]] %>%
  as.data.frame() %>%
  mutate(
    stratum = factor(stratum, levels = STRATUM_LEVELS),
    activity_binary = factor(activity_binary, levels = c("sedentary", "driving"))
  )

panel_plot <- panel_values %>%
  mutate(
    stratum = factor(stratum, levels = STRATUM_LEVELS),
    activity_binary = factor(activity_binary, levels = c("sedentary", "driving")),
    activity_label = factor(
      unname(ACTIVITY_LABELS[as.character(activity_binary)]),
      levels = unname(ACTIVITY_LABELS)
    )
  )

base_theme <- theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#EEF3F8", color = "#5E6A75", linewidth = 0.35),
    strip.text.x = element_text(face = "bold", color = "#36454F"),
    strip.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "plain", size = 11),
    plot.margin = margin(4, 4, 4, 4)
  )

# A: one row, three columns (Overall, Weekdays, Weekends)
p_dist <- ggplot(density_df, aes(x = nhr, fill = activity_binary, color = activity_binary)) +
  geom_density(alpha = 0.26, adjust = 1.0, linewidth = 0.80) +
  facet_grid(. ~ stratum) +
  scale_fill_manual(values = PAL_ACTIVITY) +
  scale_color_manual(values = PAL_ACTIVITY) +
  labs(title = "Distribution of NHR", x = "NHR [bpm]", y = "PDF") +
  base_theme

# B: occupancy above baseline
p_pct <- ggplot(panel_plot, aes(x = activity_label, y = mean_subject_pct_nhr_gt_zero, fill = activity_binary)) +
  geom_col(width = 0.68, alpha = BAR_ALPHA, color = NA) +
  geom_errorbar(
    aes(ymin = q05_subject_pct_nhr_gt_zero, ymax = q95_subject_pct_nhr_gt_zero),
    width = 0.16, linewidth = 0.70, color = COL_ERRBAR
  ) +
  facet_grid(. ~ stratum) +
  scale_fill_manual(values = PAL_ACTIVITY) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1), limits = c(0, 1.03),
    breaks = c(0, .25, .50, .75, 1.00), expand = expansion(mult = c(0, .01))
  ) +
  labs(title = "Time above baseline", x = NULL, y = "% time NHR > 0") +
  base_theme

# C: within-participant variability
sd_ymax <- max(panel_plot$q95_subject_sd_nhr, na.rm = TRUE) * 1.08
p_sd <- ggplot(panel_plot, aes(x = activity_label, y = mean_subject_sd_nhr, fill = activity_binary)) +
  geom_col(width = 0.68, alpha = BAR_ALPHA, color = NA) +
  geom_errorbar(
    aes(ymin = q05_subject_sd_nhr, ymax = q95_subject_sd_nhr),
    width = 0.16, linewidth = 0.70, color = COL_ERRBAR
  ) +
  facet_grid(. ~ stratum) +
  scale_fill_manual(values = PAL_ACTIVITY) +
  coord_cartesian(ylim = c(0, sd_ymax)) +
  labs(title = "Variability of NHR", x = NULL, y = "SD(NHR) [bpm]") +
  base_theme

# D: transition dynamics / persistence
ret_ymax <- max(panel_plot$q95_subject_returns_per_hour, na.rm = TRUE) * 1.08
p_returns <- ggplot(panel_plot, aes(x = activity_label, y = mean_subject_returns_per_hour, fill = activity_binary)) +
  geom_col(width = 0.68, alpha = BAR_ALPHA, color = NA) +
  geom_errorbar(
    aes(ymin = q05_subject_returns_per_hour, ymax = q95_subject_returns_per_hour),
    width = 0.16, linewidth = 0.70, color = COL_ERRBAR
  ) +
  facet_grid(. ~ stratum) +
  scale_fill_manual(values = PAL_ACTIVITY) +
  coord_cartesian(ylim = c(0, ret_ymax)) +
  labs(title = "Returns to baseline", x = NULL, y = "Events per observed hour") +
  base_theme

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig_expanded <- p_dist / p_pct / p_sd / p_returns +
    patchwork::plot_layout(heights = c(1.25, 1, 1, 1)) +
    patchwork::plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(
        face = "bold",
        size = 18
      ),
      plot.tag.position = c(0.008, 0.992)
    )

  safe_save_pdf(
    fig_expanded,
    file.path(fig_dir, "Figure3_Operating_Regimes.pdf"),
    w = 11.2,
    h = 11.0
  )
  safe_save_png(
    fig_expanded,
    file.path(fig_dir, "Figure3_Operating_Regimes.png"),
    w = 11.2,
    h = 11.0,
    dpi = 300
  )
} else {
  log_msg("Package patchwork is not installed; saving the four rows separately.")
  safe_save_pdf(p_dist, file.path(fig_dir, "Figure3_RowA_NHR_Distribution.pdf"), w = 11.2, h = 2.8)
  safe_save_pdf(p_pct, file.path(fig_dir, "Figure3_RowB_TimeAboveBaseline.pdf"), w = 11.2, h = 2.5)
  safe_save_pdf(p_sd, file.path(fig_dir, "Figure3_RowC_NHR_SD.pdf"), w = 11.2, h = 2.5)
  safe_save_pdf(p_returns, file.path(fig_dir, "Figure3_RowD_ReturnsToBaseline.pdf"), w = 11.2, h = 2.5)
  safe_save_png(p_dist, file.path(fig_dir, "Figure3_RowA_NHR_Distribution.png"), w = 11.2, h = 2.8, dpi = 300)
  safe_save_png(p_pct, file.path(fig_dir, "Figure3_RowB_TimeAboveBaseline.png"), w = 11.2, h = 2.5, dpi = 300)
  safe_save_png(p_sd, file.path(fig_dir, "Figure3_RowC_NHR_SD.png"), w = 11.2, h = 2.5, dpi = 300)
  safe_save_png(p_returns, file.path(fig_dir, "Figure3_RowD_ReturnsToBaseline.png"), w = 11.2, h = 2.5, dpi = 300)
}

# -------------------------------------------------------------------------

log_msg("Saved figures under: ", fig_dir)
log_msg("Saved manuscript Table 3 under: ", table_dir)
log_msg("Session information:")
capture.output(sessionInfo(), file = log_file, append = TRUE)
log_msg("DONE.")
