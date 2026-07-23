# ============================================================
# 02_figure2_representative_trace.R
#
# PURPOSE
#   Generate Figure 2 for the revised npj Digital Public Health manuscript:
#
#     Wearable sensing reveals the structure of cardiac
#     activation associated with everyday driving
#
#   Figure 2 presents a representative participant-level heart-rate
#   trace across study days 1-7. The figure shows:
#
#     - observed raw heart rate, HR_raw,itd, colored by behavioral context;
#     - the available participant-day baseline, HR_base,id, overlaid in red.
#
#   The displayed participant-day baseline is used for visualization only.
#   Subsequent baseline-referenced analyses in the manuscript use the
#   participant-level baseline, HR_base,i, defined as the mean of each
#   participant's available daily baseline values.
#
# INPUT
#   Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
# DEFAULT SETTINGS
#   Subject:    P62
#   Resolution: 60 s
#
#   The defaults can be changed without editing the script by setting:
#
#     NUBI_SUBJECT_ID=P62
#     NUBI_RES_SECONDS=60
#
# REQUIRED ACTIVITY LABELS
#   activity3 == "driving"
#   activity3 == "non_driving_sedentary"
#   activity3 == "non_driving_physical_activity"
#
# STUDY-DAY AXIS
#   day_num follows the Monday-Sunday study-week convention described in
#   the Methods: day 1 is Monday and day 7 is Sunday. Each study day is
#   placed in a fixed 24-h slot, yielding a 0-168 h display.
#
# MAJOR OUTPUTS
#   Results/paper_figs/<timestamp>_<resolution>sec_figure2_<subject>/
#
#     Figures/Figure2_<subject>.pdf
#     Figures/Figure2_<subject>.png
#     Diagnostics/activity3_levels_detected.csv
#     Diagnostics/counts_<subject>.csv
#     Diagnostics/subject_summary_<subject>.csv
#     Diagnostics/baseline_study_day_<subject>.csv
#     Diagnostics/missingness_after_fill_<subject>.csv
#     Diagnostics/duplicate_pid_time_<subject>.csv
#     run_log.txt
#
# REPOSITORY SCOPE
#   The public repository starts from the curated analysis-ready dataset.
#   It does not reconstruct the dataset from raw wearable, smartphone,
#   vehicle, or ground-truth streams.
#
# PRIVACY NOTE
#   Direct GPS coordinates are not required for this figure.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(ggplot2)
})

options(warn = 1)

# ----------------------------
# User-configurable settings
# ----------------------------
LOCAL_TZ <- "America/Chicago"

SUBJECT_ID <- trimws(Sys.getenv("NUBI_SUBJECT_ID", unset = "P62"))
if (!nzchar(SUBJECT_ID)) stop("NUBI_SUBJECT_ID cannot be empty.")

RES_SECONDS <- suppressWarnings(
  as.integer(Sys.getenv("NUBI_RES_SECONDS", unset = "60"))
)
if (is.na(RES_SECONDS) || !RES_SECONDS %in% c(10L, 30L, 60L)) {
  stop("NUBI_RES_SECONDS must be one of 10, 30, or 60.")
}

KEEP_DAYS <- 7L
X_END_HRS <- 24 * KEEP_DAYS
X_BREAK_BY_HRS <- 6L
# Small right-side plotting allowance prevents the final "24:00 / 7" tick label
# from being clipped in PDF and PNG output.
X_RIGHT_PAD_HRS <- 1.5
GAP_MULT <- 2L
GAP_THRESHOLD_SECS <- as.integer(GAP_MULT * RES_SECONDS)

PDF_W <- 11
PDF_H <- 4.8
PNG_DPI <- 300
BASE_FONT <- 12

# Validated activity3 labels in the curated dataset.
ACT3_DRIVING <- "driving"
ACT3_ND_SED  <- "non_driving_sedentary"
ACT3_ND_PA   <- "non_driving_physical_activity"

# Internal plotting states use manuscript terminology.
STATE_LEVELS <- c("DRIVING", "NONDRIVING_SEDENTARY", "PHYSICAL_ACTIVITY")

PAL_STATE <- c(
  "DRIVING" = "orange",
  "NONDRIVING_SEDENTARY" = "black",
  "PHYSICAL_ACTIVITY" = "springgreen3",
  "Participant-day baseline" = "red"
)

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

# ----------------------------
# Input and output paths
# ----------------------------
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
  "Results",
  "paper_figs",
  paste0(stamp, "_", RES_SECONDS, "sec_figure2_", SUBJECT_ID)
)
fig_dir <- file.path(out_dir, "Figures")
diag_dir <- file.path(out_dir, "Diagnostics")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(fig_dir, paste0("Figure2_", SUBJECT_ID, ".pdf"))
png_path <- file.path(fig_dir, paste0("Figure2_", SUBJECT_ID, ".png"))
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

# ============================================================
# Helpers
# ============================================================
norm_activity <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x <- gsub("\\s+", "_", x)
  gsub("-", "_", x)
}

parse_time_local <- function(x, tz = LOCAL_TZ) {
  if (inherits(x, "POSIXt")) return(force_tz(x, tzone = tz))

  x2 <- trimws(as.character(x))
  x2 <- gsub("Z$", "", x2, ignore.case = TRUE)

  tt <- suppressWarnings(ymd_hms(x2, tz = tz, quiet = TRUE))
  if (!all(is.na(tt))) return(tt)

  suppressWarnings(parse_date_time(
    x2,
    orders = c(
      "ymd HMS", "ymd HM", "ymdT HMS", "ymdT HM",
      "mdy HMS", "mdy HM", "dmy HMS", "dmy HM",
      "Ymd HMS", "Ymd HM", "YmdT HMS", "YmdT HM",
      "ymd HMSOS", "mdy HMSOS", "dmy HMSOS"
    ),
    tz = tz
  ))
}

parse_day_num <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x[x %in% c("", "na", "n/a", "null", "unknown")] <- NA_character_
  out <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  out[!out %in% seq_len(KEEP_DAYS)] <- NA_integer_
  out
}

map_activity_state <- function(x) {
  a <- norm_activity(x)
  out <- fifelse(
    a == ACT3_DRIVING,
    "DRIVING",
    fifelse(
      a == ACT3_ND_SED,
      "NONDRIVING_SEDENTARY",
      fifelse(a == ACT3_ND_PA, "PHYSICAL_ACTIVITY", NA_character_)
    )
  )
  factor(out, levels = STATE_LEVELS)
}

mode_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

pad_range <- function(r, frac = 0.06) {
  if (length(r) != 2L || any(!is.finite(r))) return(c(0, 1))
  span <- diff(r)
  if (!is.finite(span) || span <= 0) return(r + c(-1, 1))
  r + c(-1, 1) * span * frac
}

clock_labels_7day <- function(x) {
  out <- character(length(x))
  last_x <- X_END_HRS - RES_SECONDS / 3600
  near_end <- abs(x - last_x) < 1e-9
  out[near_end] <- "24:00\n7"

  xx <- x[!near_end]
  day_idx <- floor(xx / 24) + 1L
  h <- ((xx %% 24) + 24) %% 24
  hh <- floor(h)
  mm <- round((h - hh) * 60)
  mm[mm == 60] <- 0
  hh[hh == 24] <- 0

  out[!near_end] <- paste0(sprintf("%02d:%02d", hh, mm), "\n", day_idx)
  out
}

fill_missing_grid_one_subject <- function(d, step_secs) {
  setorder(d, day_num_int, sec_of_day)

  out_list <- vector("list", KEEP_DAYS)

  for (day_index in seq_len(KEEP_DAYS)) {
    d_day <- d[day_num_int == day_index]

    full <- data.table(
      p_id = d$p_id[1],
      day_num_int = day_index,
      sec_of_day = seq(0, 24 * 3600 - step_secs, by = step_secs)
    )

    if (nrow(d_day) > 0L) {
      observed <- copy(d_day[, .(
        p_id,
        day_num_int,
        sec_of_day,
        dt_time,
        raw_hr,
        baseline_hr,
        state
      )])

      out <- merge(
        full,
        observed,
        by = c("p_id", "day_num_int", "sec_of_day"),
        all.x = TRUE,
        sort = TRUE
      )
    } else {
      out <- copy(full)
      out[, `:=`(
        dt_time = as.POSIXct(NA, tz = LOCAL_TZ),
        raw_hr = NA_real_,
        baseline_hr = NA_real_,
        state = factor(NA_character_, levels = STATE_LEVELS)
      )]
    }

    # Display the participant-day baseline only where an observed HR value exists,
    # matching the manuscript figure rather than drawing across unobserved periods.
    out[!is.finite(raw_hr), `:=`(
      baseline_hr = NA_real_,
      state = factor(NA_character_, levels = STATE_LEVELS)
    )]

    out_list[[day_index]] <- out
  }

  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

make_subject_plot <- function(dsub) {
  setorder(dsub, x_hr)

  gap_threshold_hours <- GAP_THRESHOLD_SECS / 3600
  dsub[, dx_hr := x_hr - shift(x_hr)]

  # Raw-HR segments break at temporal gaps, missingness boundaries, day
  # boundaries, and behavioral-context transitions.
  dsub[, new_raw_segment :=
         is.na(dx_hr) |
         dx_hr < 0 |
         dx_hr > gap_threshold_hours |
         day_num_int != shift(day_num_int) |
         state != shift(state) |
         is.na(raw_hr) != shift(is.na(raw_hr))]
  dsub[, raw_segment_id := cumsum(fifelse(is.na(new_raw_segment), TRUE, new_raw_segment))]

  # Baseline segments break only at temporal gaps, missingness boundaries,
  # and day boundaries. This preserves a continuous daily baseline across
  # changes in behavioral context.
  dsub[, new_baseline_segment :=
         is.na(dx_hr) |
         dx_hr < 0 |
         dx_hr > gap_threshold_hours |
         day_num_int != shift(day_num_int) |
         is.na(baseline_hr) != shift(is.na(baseline_hr))]
  dsub[, baseline_segment_id := cumsum(
    fifelse(is.na(new_baseline_segment), TRUE, new_baseline_segment)
  )]

  y_values <- c(dsub$raw_hr, dsub$baseline_hr)
  y_values <- y_values[is.finite(y_values)]
  y_limits <- if (length(y_values)) pad_range(range(y_values)) else c(0, 1)

  ggplot(as.data.frame(dsub), aes(x = x_hr)) +
    geom_line(
      aes(y = raw_hr, color = state, group = raw_segment_id),
      linewidth = 0.38,
      na.rm = TRUE
    ) +
    geom_line(
      aes(
        y = baseline_hr,
        color = "Participant-day baseline",
        group = baseline_segment_id
      ),
      linewidth = 0.48,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = PAL_STATE,
      breaks = c(
        "DRIVING",
        "NONDRIVING_SEDENTARY",
        "PHYSICAL_ACTIVITY",
        "Participant-day baseline"
      ),
      labels = c(
        "Driving",
        "Non-driving sedentary",
        "Physical activity",
        "Participant-day baseline"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(
        seq(0, X_END_HRS - X_BREAK_BY_HRS, by = X_BREAK_BY_HRS),
        X_END_HRS - RES_SECONDS / 3600
      ),
      labels = clock_labels_7day,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(
        0,
        X_END_HRS - RES_SECONDS / 3600 + X_RIGHT_PAD_HRS
      ),
      ylim = y_limits,
      clip = "off"
    ) +
    labs(
      title = paste0("Participant ", SUBJECT_ID, ": Raw HR across study days 1–7"),
      subtitle = paste0("Resolution: ", RES_SECONDS, " sec"),
      x = NULL,
      y = "HR [bpm]",
      color = NULL
    ) +
    theme_minimal(base_size = BASE_FONT) +
    theme(
      plot.title = element_text(face = "bold", size = BASE_FONT + 1),
      plot.subtitle = element_text(size = BASE_FONT - 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
    )
}

# ============================================================
# Main analysis
# ============================================================
main <- function() {
  log_msg("SCRIPT VERSION: GITHUB-READY MANUSCRIPT-SYNCHRONIZED FINAL 2026-07-23")
  log_msg("Script: 02_figure2_representative_trace.R")
  log_msg("Subject: ", SUBJECT_ID)
  log_msg("Resolution: ", RES_SECONDS, " sec")
  log_msg("Input: ", in_path)
  log_msg("Output directory: ", out_dir)
  log_msg("Reading dataset...")

  dt <- fread(in_path, showProgress = TRUE)
  log_msg("Rows read: ", format(nrow(dt), big.mark = ","))
  log_msg("Columns read: ", ncol(dt))

  required_columns <- c("p_id", "time", "raw_hr", "bl_hr", "activity3", "day_num")
  missing_columns <- setdiff(required_columns, names(dt))
  if (length(missing_columns)) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  dt[, p_id := as.character(p_id)]
  dt[, time := as.character(time)]
  dt[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]

  n_bad_time <- dt[is.na(dt_time), .N]
  if (n_bad_time > 0L) log_msg("Rows with unparseable time removed: ", n_bad_time)
  dt <- dt[!is.na(dt_time)]

  dt[, `:=`(
    raw_hr = suppressWarnings(as.numeric(raw_hr)),
    baseline_hr = suppressWarnings(as.numeric(bl_hr)),
    activity3_norm = norm_activity(activity3),
    day_num_int = parse_day_num(day_num)
  )]

  invalid_day_rows <- dt[is.na(day_num_int), .N]
  if (invalid_day_rows > 0L) {
    log_msg("Rows with invalid study-day labels removed: ", invalid_day_rows)
  }
  dt <- dt[!is.na(day_num_int)]

  activity_levels <- dt[, .N, by = activity3_norm][order(-N)]
  fwrite(activity_levels, file.path(diag_dir, "activity3_levels_detected.csv"))

  dt[, state := map_activity_state(activity3_norm)]
  if (anyNA(dt$state)) {
    unknown <- sort(unique(dt[is.na(state), activity3_norm]))
    stop("Unrecognized activity3 label(s): ", paste(unknown, collapse = ", "))
  }

  # Keep only the requested participant after global label validation.
  dt <- dt[p_id == SUBJECT_ID]
  if (!nrow(dt)) stop("No rows found for subject: ", SUBJECT_ID)

  duplicate_rows <- dt[, .N, by = .(p_id, dt_time)][N > 1L]
  fwrite(
    duplicate_rows,
    file.path(diag_dir, paste0("duplicate_pid_time_", SUBJECT_ID, ".csv"))
  )
  log_msg("Duplicate (participant, timestamp) groups: ", nrow(duplicate_rows))

  if (nrow(duplicate_rows) > 0L) {
    state_levels <- STATE_LEVELS
    dt <- dt[, .(
      raw_hr = if (all(!is.finite(raw_hr))) NA_real_ else mean(raw_hr, na.rm = TRUE),
      baseline_hr = if (all(!is.finite(baseline_hr))) NA_real_ else mean(baseline_hr, na.rm = TRUE),
      state = factor(mode_chr(state), levels = state_levels),
      day_num_int = suppressWarnings(as.integer(mode_chr(day_num_int)))
    ), by = .(p_id, dt_time)]
  } else {
    dt <- dt[, .(p_id, dt_time, raw_hr, baseline_hr, state, day_num_int)]
  }

  dt[, sec_of_day :=
       hour(dt_time) * 3600L +
       minute(dt_time) * 60L +
       floor(second(dt_time))]

  # Require alignment to the selected resolution. A mismatch indicates that
  # the chosen curated file and the requested resolution are inconsistent.
  off_grid <- dt[sec_of_day %% RES_SECONDS != 0L, .N]
  if (off_grid > 0L) {
    log_msg("Observed rows not aligned to the ", RES_SECONDS, "-s grid: ", off_grid)
  }

  counts <- dt[, .N, by = .(p_id, state)]
  counts_wide <- dcast(counts, p_id ~ state, value.var = "N", fill = 0)
  for (context_name in STATE_LEVELS) {
    if (!context_name %in% names(counts_wide)) {
      counts_wide[, (context_name) := 0L]
    }
  }
  counts_wide[, TOTAL :=
                DRIVING + NONDRIVING_SEDENTARY + PHYSICAL_ACTIVITY]
  setcolorder(counts_wide, c("p_id", STATE_LEVELS, "TOTAL"))
  fwrite(
    counts_wide,
    file.path(diag_dir, paste0("counts_", SUBJECT_ID, ".csv"))
  )

  subject_summary <- dt[, .(
    n_rows = .N,
    n_study_days_present = uniqueN(day_num_int),
    days_present = paste(sort(unique(day_num_int)), collapse = ","),
    time_min = min(dt_time),
    time_max = max(dt_time),
    n_raw_hr_missing = sum(!is.finite(raw_hr)),
    pct_raw_hr_missing = 100 * mean(!is.finite(raw_hr)),
    n_baseline_missing = sum(!is.finite(baseline_hr)),
    pct_baseline_missing = 100 * mean(!is.finite(baseline_hr))
  ), by = p_id]
  fwrite(
    subject_summary,
    file.path(diag_dir, paste0("subject_summary_", SUBJECT_ID, ".csv"))
  )

  baseline_diagnostic <- dt[, .(
    n_rows_observed = .N,
    any_raw_hr = any(is.finite(raw_hr)),
    any_baseline = any(is.finite(baseline_hr)),
    n_nonmissing_baseline_rows = sum(is.finite(baseline_hr)),
    n_unique_baseline_values = uniqueN(baseline_hr[is.finite(baseline_hr)]),
    baseline_values_seen = paste(
      sort(unique(baseline_hr[is.finite(baseline_hr)])),
      collapse = ";"
    )
  ), by = .(p_id, day_num_int)][order(day_num_int)]
  fwrite(
    baseline_diagnostic,
    file.path(diag_dir, paste0("baseline_study_day_", SUBJECT_ID, ".csv"))
  )

  multiple_baseline_days <- baseline_diagnostic[n_unique_baseline_values > 1L]
  if (nrow(multiple_baseline_days) > 0L) {
    stop(
      "More than one participant-day baseline value was detected for subject ",
      SUBJECT_ID, " on study day(s): ",
      paste(multiple_baseline_days$day_num_int, collapse = ", ")
    )
  }

  log_msg("Rows kept for subject before fixed-grid expansion: ", nrow(dt))
  log_msg(
    "Study days represented: ",
    paste(sort(unique(dt$day_num_int)), collapse = ", ")
  )

  dt_filled <- fill_missing_grid_one_subject(dt, step_secs = RES_SECONDS)
  setorder(dt_filled, p_id, day_num_int, sec_of_day)
  dt_filled[, x_hr :=
              24 * (day_num_int - 1L) + sec_of_day / 3600]

  fill_diagnostic <- dt_filled[, .(
    n_total_grid_points = .N,
    n_observed_raw_hr = sum(is.finite(raw_hr)),
    n_missing_raw_hr = sum(!is.finite(raw_hr)),
    pct_missing_raw_hr = 100 * mean(!is.finite(raw_hr))
  ), by = p_id]
  fwrite(
    fill_diagnostic,
    file.path(diag_dir, paste0("missingness_after_fill_", SUBJECT_ID, ".csv"))
  )

  expected_grid_rows <- KEEP_DAYS * 24L * 3600L / RES_SECONDS
  if (nrow(dt_filled) != expected_grid_rows) {
    stop(
      "Fixed-grid construction failed: expected ", expected_grid_rows,
      " rows but obtained ", nrow(dt_filled), "."
    )
  }

  log_msg("Rows after fixed-grid expansion: ", nrow(dt_filled))

  figure2 <- make_subject_plot(dt_filled)

  safe_save_pdf(figure2, pdf_path, w = PDF_W, h = PDF_H)
  ggsave(
    filename = png_path,
    plot = figure2,
    width = PDF_W,
    height = PDF_H,
    dpi = PNG_DPI
  )

  log_msg("Wrote Figure 2 PDF: ", pdf_path)
  log_msg("Wrote Figure 2 PNG: ", png_path)
  log_msg("Session information:")
  capture.output(sessionInfo(), file = log_file, append = TRUE)
  log_msg("DONE.")
}

main()
