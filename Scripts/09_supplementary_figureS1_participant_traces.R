# ============================================================
# 09_supplementary_figureS1_participant_traces_v1.0.R
#
# PURPOSE
#   Reproduce Supplementary Figure S1 for the npj Digital Public Health
#   manuscript:
#
#     Wearable sensing reveals the structure of cardiac activation
#     associated with everyday driving
#
#   The figure visualizes participant-level heart-rate time series
#   across the 7 study-day slots, with raw HR colored by behavioral
#   context and baseline HR overlaid.
#
# DATA SOURCE
#   Final clean MASTER dataset selected by RES_SECONDS:
#
#     Data/NUBI_Data_10sec_Level_MASTER_CLEAN.csv
#     Data/NUBI_Data_30sec_Level_MASTER_CLEAN.csv
#     Data/NUBI_Data_60sec_Level_MASTER_CLEAN.csv
#
#   RES_SECONDS defaults to 60, which reproduces the manuscript figure.
#   Advanced users may override it with the NUBI_RES_SECONDS environment
#   variable (10, 30, or 60), provided the matching MASTER file is present.
#
# IMPORTANT COLUMN NOTE
#   The MASTER dataset contains both:
#
#     activity  = two-class label: driving / not_driving
#     activity3 = three-class label:
#                 driving
#                 non_driving_sedentary
#                 non_driving_physical_activity
#
#   This figure MUST use activity3, not activity.
#
# PANEL DEFINITION
#   One subject per panel, with multiple subjects per page.
#
#   Raw HR is colored by behavioral context derived from activity3:
#
#     driving                         -> orange
#     non_driving_sedentary           -> black
#     non_driving_physical_activity   -> green
#
#   Baseline HR is shown in red.
#
# X-AXIS DEFINITION
#   The x-axis uses study-day placement over a fixed 7-slot layout:
#
#     day1 = Monday slot, 0-24 h
#     day2 = Tuesday slot, 24-48 h
#     ...
#     day7 = Sunday slot, 144-168 h
#
#   day_num is treated as the primary study-day key. Day1..Day7 are
#   weekday-coded study-day labels and are not assumed to be
#   consecutive calendar dates.
#
# TIME HANDLING
#   Timestamps in `time` are parsed as local clock time in:
#
#     America/Chicago
#
#   Missing grid points are inserted within each study-day slot at
#   the selected resolution so gaps appear as line breaks.
#
# REQUIRED COLUMNS
#   After name canonicalization, the input dataset must contain:
#
#     p_id
#     time
#     raw_hr
#     bl_hr
#     activity3
#     day_num
#
# MAJOR OUTPUTS
#   The script writes the supplementary figure and diagnostics under:
#
#     Results/paper_figs/<timestamp>_<RES>sec_supplementary_figureS1_participant_traces/
#
#   Main figure:
#
#     Figures/Supplementary_Figure_S1.pdf
#
#   Diagnostics:
#
#     Diagnostics/activity3_counts.csv
#     Diagnostics/mapped_state_counts.csv
#     Diagnostics/unrecognized_activity3_labels.csv, only if needed
#     Diagnostics/subject_summary_before_fill.csv
#     Diagnostics/duplicate_pid_time_rows.csv
#     Diagnostics/subject_summary_study_day_slots.csv
#     Diagnostics/baseline_subject_study_day_level.csv
#     Diagnostics/baseline_subject_study_day_summary.csv
#     Diagnostics/baseline_study_day_overall_summary.csv
#     Diagnostics/baseline_subject_calendar_day_level.csv
#     Diagnostics/baseline_subject_calendar_day_summary.csv
#     Diagnostics/baseline_calendar_day_overall_summary.csv
#     Diagnostics/missingness_after_fill.csv
#     run_log.txt
#     session_info.txt
#
# REPOSITORY SCOPE
#   This public repository starts from the final clean MASTER
#   dataset. It does not rebuild upstream wearable, smartphone,
#   vehicle, or ground-truth streams.
#
# PRIVACY NOTE
#   This script uses curated timestamped features and activity
#   labels. It does not require direct GPS coordinates or raw
#   location traces.
# ============================================================

required_packages <- c("data.table", "lubridate", "stringr", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

options(warn = 1)

# ----------------------------
# USER SETTINGS
# ----------------------------
SCRIPT_VERSION <- "1.0"
LOCAL_TZ <- "America/Chicago"

subjects_per_page <- 3
x_break_by <- 6
use_global_y <- TRUE

KEEP_DAYS <- 7L
X_END_HRS <- 24 * KEEP_DAYS

# Colors
COL_DRIVE <- "orange"
COL_SED   <- "black"
COL_ACT   <- "springgreen3"
COL_BASE  <- "red"

STATE_LEVELS <- c(
  "DRIVING",
  "NONDRIVING_SEDENTARY",
  "NONDRIVING_PHYSICAL_ACTIVITY"
)

PAL_STATE <- c(
  "DRIVING" = COL_DRIVE,
  "NONDRIVING_SEDENTARY" = COL_SED,
  "NONDRIVING_PHYSICAL_ACTIVITY" = COL_ACT
)

# ----------------------------
# Resolution and project root
# ----------------------------
# The manuscript figure uses 60-sec data. An environment-variable override is
# retained for optional exploratory rendering at 10 or 30 sec without making
# the public reproduction workflow interactive.
RES_SECONDS <- suppressWarnings(as.integer(
  Sys.getenv("NUBI_RES_SECONDS", unset = "60")
))
if (length(RES_SECONDS) != 1L || is.na(RES_SECONDS) ||
    !(RES_SECONDS %in% c(10L, 30L, 60L))) {
  stop("NUBI_RES_SECONDS must be one of 10, 30, or 60.")
}

file_map <- c(
  "10" = "NUBI_Data_10sec_Level_MASTER_CLEAN.csv",
  "30" = "NUBI_Data_30sec_Level_MASTER_CLEAN.csv",
  "60" = "NUBI_Data_60sec_Level_MASTER_CLEAN.csv"
)

in_file <- unname(file_map[as.character(RES_SECONDS)])

PROJECT_ROOT <- Sys.getenv("NUBI_PROJECT_ROOT", unset = NA_character_)
if (is.na(PROJECT_ROOT) || PROJECT_ROOT == "") {
  candidate_roots <- unique(normalizePath(
    c(getwd(), file.path(getwd(), "..")),
    mustWork = FALSE
  ))
  has_data_file <- file.exists(file.path(candidate_roots, "Data", in_file))
  if (!any(has_data_file)) {
    stop(
      "Could not locate Data/", in_file, ". ",
      "Run this script from the project root or Scripts/ folder, or set ",
      "Sys.setenv(NUBI_PROJECT_ROOT = '/path/to/project')."
    )
  }
  PROJECT_ROOT <- candidate_roots[which(has_data_file)[1]]
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, mustWork = TRUE)

in_path <- file.path(PROJECT_ROOT, "Data", in_file)
if (!file.exists(in_path)) {
  stop("Chosen input file does not exist: ", in_path)
}

# Break lines across time gaps > 2 expected bins
gap_threshold_secs <- 2L * RES_SECONDS

# ----------------------------
# Output folders
# ----------------------------
stamp   <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(
  PROJECT_ROOT, "Results", "paper_figs",
  paste0(stamp, "_", RES_SECONDS, "sec_supplementary_figureS1_participant_traces")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

diag_dir <- file.path(out_dir, "Diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf   <- file.path(fig_dir, "Supplementary_Figure_S1.pdf")
log_file  <- file.path(out_dir, "run_log.txt")

# ----------------------------
# Helpers
# ----------------------------
snakeify <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

canonicalize_names <- function(dt) {
  stopifnot(is.data.table(dt))
  
  old <- names(dt)
  sn  <- snakeify(old)
  
  # Important:
  # The dataset contains both:
  #   activity  = driving / not_driving
  #   activity3 = driving / non_driving_sedentary / non_driving_physical_activity
  #
  # For this figure, preserve the true activity3 column.
  # Do NOT canonicalize activity -> activity3 when activity3 is present.
  
  canon <- c(
    p_id = "p_id",
    pid = "p_id",
    participant_id = "p_id",
    participantid = "p_id",
    
    time = "time",
    timestamp = "time",
    datetime = "time",
    date_time = "time",
    
    raw_hr = "raw_hr",
    hr = "raw_hr",
    
    bl_hr = "bl_hr",
    baseline_hr = "bl_hr",
    hr_bl = "bl_hr",
    hrbl = "bl_hr",
    
    activity3 = "activity3",
    activity_3 = "activity3",
    
    day_num = "day_num",
    daynum = "day_num",
    days = "day_num"
  )
  
  new <- sn
  hit <- sn %in% names(canon)
  new[hit] <- unname(canon[sn[hit]])
  
  # Fallback only for older files that lack activity3 entirely.
  # This should NOT trigger for the current MASTER file.
  has_activity3 <- any(new == "activity3")
  if (!has_activity3 && any(sn %in% c("activity", "activity_type", "context"))) {
    idx <- which(sn %in% c("activity", "activity_type", "context"))[1]
    new[idx] <- "activity3"
  }
  
  if (!identical(old, new)) {
    new <- make.unique(new, sep = "_")
    setnames(dt, old, new)
  }
  
  dt
}

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

open_pdf <- function(path, width, height) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(path, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::pdf(path, width = width, height = height, onefile = TRUE, useDingbats = FALSE)
  }
}

pad_range <- function(r, frac = 0.06) {
  if (length(r) != 2 || any(!is.finite(r))) return(c(0, 1))
  span <- diff(r)
  if (!is.finite(span) || span <= 0) return(r + c(-1, 1))
  r + c(-1, 1) * span * frac
}

clock_labels_7day <- function(x) {
  out <- character(length(x))
  
  near_end <- abs(x - (X_END_HRS - RES_SECONDS / 3600)) < 1e-9
  out[near_end] <- "24:00\n7"
  
  xx <- x[!near_end]
  day_idx <- floor(xx / 24) + 1
  h <- ((xx %% 24) + 24) %% 24
  hh <- floor(h)
  mm <- round((h - hh) * 60)
  mm[mm == 60] <- 0
  hh[hh == 24] <- 0
  hhmm <- sprintf("%02d:%02d", hh, mm)
  out[!near_end] <- paste0(hhmm, "\n", day_idx)
  
  out
}

norm_chr <- function(x) {
  trimws(tolower(as.character(x)))
}

activity_key <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

parse_time_local <- function(x, tz = LOCAL_TZ) {
  if (inherits(x, "POSIXt")) return(force_tz(x, tzone = tz))
  
  x <- as.character(x)
  x2 <- sub("Z$", "", x)
  
  tt <- suppressWarnings(ymd_hms(x2, tz = tz, quiet = TRUE))
  if (sum(!is.na(tt)) > 0) return(tt)
  
  suppressWarnings(parse_date_time(
    x2,
    orders = c(
      "ymd HMS", "ymd HM",
      "mdy HMS", "mdy HM",
      "dmy HMS", "dmy HM",
      "ymd HMSOS", "mdy HMSOS", "dmy HMSOS"
    ),
    tz = tz
  ))
}

pid_sort_key <- function(pid_chr) {
  pid_chr <- as.character(pid_chr)
  dig <- suppressWarnings(as.integer(str_extract(pid_chr, "\\d+")))
  dig[is.na(dig)] <- 1e9L
  dig
}

mode_chr <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

state_from_activity3 <- function(activity3_vec) {
  k <- activity_key(activity3_vec)
  
  out <- rep(NA_character_, length(k))
  
  out[k == "driving"] <- "DRIVING"
  out[k == "nondrivingsedentary"] <- "NONDRIVING_SEDENTARY"
  out[k == "nondrivingphysicalactivity"] <- "NONDRIVING_PHYSICAL_ACTIVITY"
  
  factor(out, levels = STATE_LEVELS)
}

parse_day_num <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "n/a", "null")] <- NA_character_
  out <- suppressWarnings(as.integer(str_extract(x, "\\d+")))
  out[!(out %in% 1:7)] <- NA_integer_
  out
}

write_diag_table <- function(dt, path) {
  fwrite(dt, path)
}

make_page_header <- function(pg, n_pg) {
  leg <- data.frame(
    x = c(0.02, 0.30, 0.63),
    y = c(0.25, 0.25, 0.25),
    col = c(COL_DRIVE, COL_SED, COL_ACT),
    lab = c("Driving", "Non-driving sedentary", "Non-driving physical activity")
  )
  
  ggplot() +
    theme_void(base_size = 11) +
    annotate(
      "text",
      x = 0,
      y = 0.98,
      hjust = 0,
      vjust = 1,
      size = 5,
      fontface = "bold",
      label = "Supplementary Figure"
    ) +
    annotate(
      "text",
      x = 1,
      y = 0.98,
      hjust = 1,
      vjust = 1,
      size = 3.8,
      label = paste0("Resolution: ", RES_SECONDS, " sec | Page ", pg, " of ", n_pg)
    ) +
    annotate(
      "text",
      x = 0,
      y = 0.83,
      hjust = 0,
      vjust = 1,
      size = 3.6,
      label = paste0(
        "Raw HR colored by behavioral context (", RES_SECONDS,
        "-sec). Baseline in red. X-axis uses study-day slots day1..day7; missing samples are blank."
      )
    ) +
    geom_point(
      data = leg,
      aes(x = x, y = y),
      color = leg$col,
      size = 3.2,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = leg,
      aes(x = x + 0.03, y = y, label = lab),
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      inherit.aes = FALSE
    ) +
    geom_segment(
      aes(x = 0.86, xend = 0.92, y = 0.25, yend = 0.25),
      color = COL_BASE,
      linewidth = 1.1,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = 0.93, y = 0.25, label = "Baseline"),
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      inherit.aes = FALSE
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")
}

make_subject_plot <- function(dsub, subject_label, global_y = NULL) {
  setorder(dsub, x_hr, dt_time)
  
  dsub[, dx_hr := x_hr - shift(x_hr, 1)]
  
  gap_threshold_hours <- gap_threshold_secs / 3600
  
  dsub[, new_seg := is.na(dx_hr) |
         (dx_hr < 0) |
         (dx_hr > gap_threshold_hours) |
         (state != shift(state, 1)) |
         (is.na(raw_hr) != shift(is.na(raw_hr), 1))
  ]
  dsub[, seg_id := cumsum(new_seg)]
  
  if (!is.null(global_y) && all(is.finite(global_y))) {
    y_lim <- global_y
  } else {
    y_vals <- c(dsub$raw_hr, dsub$baseline_hr)
    y_vals <- y_vals[is.finite(y_vals)]
    y_lim <- if (length(y_vals) > 0) pad_range(range(y_vals)) else c(0, 1)
  }
  
  ggplot(as.data.frame(dsub), aes(x = x_hr)) +
    geom_line(
      aes(y = raw_hr, color = state, group = seg_id),
      linewidth = 0.40,
      lineend = "round",
      na.rm = TRUE
    ) +
    geom_line(
      aes(y = baseline_hr, group = seg_id),
      color = COL_BASE,
      linewidth = 0.60,
      lineend = "round",
      na.rm = TRUE
    ) +
    scale_color_manual(values = PAL_STATE, drop = FALSE) +
    scale_x_continuous(
      breaks = c(
        seq(0, X_END_HRS - x_break_by, by = x_break_by),
        X_END_HRS - RES_SECONDS / 3600
      ),
      labels = clock_labels_7day,
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(0, X_END_HRS - RES_SECONDS / 3600),
      ylim = y_lim
    ) +
    labs(title = subject_label, x = NULL, y = "HR [bpm]") +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 6, 2, 6)
    )
}

fill_missing_grid_one_subject <- function(d, step_secs) {
  setorder(d, day_num_int, dt_time)
  
  slots <- 1:KEEP_DAYS
  out_list <- vector("list", length(slots))
  
  for (ii in seq_along(slots)) {
    dd <- slots[ii]
    d_day <- d[day_num_int == dd]
    
    full <- data.table(
      sec_of_day = seq(0, 24 * 3600 - step_secs, by = step_secs)
    )
    full[, `:=`(
      p_id = d$p_id[1],
      day_num_int = dd
    )]
    
    if (nrow(d_day) > 0) {
      tmp <- copy(d_day[, .(
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
        tmp,
        by = c("p_id", "day_num_int", "sec_of_day"),
        all.x = TRUE,
        sort = TRUE
      )
    } else {
      out <- copy(full)
      out[, `:=`(
        dt_time = as.POSIXct(NA),
        raw_hr = NA_real_,
        baseline_hr = NA_real_,
        state = factor(NA_character_, levels = STATE_LEVELS)
      )]
    }
    
    # Blank baseline and state at inserted grid points so missing data are
    # visually blank rather than showing a red baseline line across gaps.
    out[!is.finite(raw_hr), `:=`(
      baseline_hr = NA_real_,
      state = factor(NA_character_, levels = STATE_LEVELS)
    )]
    
    out_list[[ii]] <- out
  }
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}

main <- function() {
  log_msg("Supplementary Figure S1 generation start")
  log_msg("Script version: ", SCRIPT_VERSION)
  log_msg("Chosen resolution: ", RES_SECONDS, " sec")
  log_msg("Input: ", in_path)
  log_msg("Output dir: ", out_dir)
  log_msg("Writing PDF: ", out_pdf)
  
  dt <- fread(in_path, showProgress = TRUE)
  dt <- canonicalize_names(dt)
  
  log_msg("Rows read: ", format(nrow(dt), big.mark = ","))
  log_msg("Columns read: ", ncol(dt))
  
  need <- c("p_id", "time", "raw_hr", "bl_hr", "activity3", "day_num")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns after canonicalization: ", paste(miss, collapse = ", "))
  }
  
  dt[, p_id := as.character(p_id)]
  dt[, time := as.character(time)]
  dt[, dt_time := parse_time_local(time, tz = LOCAL_TZ)]
  
  n_bad_time <- dt[is.na(dt_time), .N]
  if (n_bad_time > 0) {
    log_msg("Rows with unparseable time removed: ", n_bad_time)
  }
  dt <- dt[!is.na(dt_time)]
  
  dt[, raw_hr := suppressWarnings(as.numeric(raw_hr))]
  dt[, baseline_hr := suppressWarnings(as.numeric(bl_hr))]
  dt[, activity3_norm := norm_chr(activity3)]
  dt[, activity3_key := activity_key(activity3)]
  
  # Confirm we are using the true three-class activity3 column.
  state_counts <- dt[, .N, by = .(activity3_norm, activity3_key)][order(-N)]
  write_diag_table(state_counts, file.path(diag_dir, "activity3_counts.csv"))
  print(state_counts)
  
  dt[, state := state_from_activity3(activity3_norm)]
  
  bad_activity <- dt[is.na(state), .N]
  if (bad_activity > 0) {
    bad_activity_labels <- dt[is.na(state), .N, by = .(activity3_norm, activity3_key)][order(-N)]
    write_diag_table(
      bad_activity_labels,
      file.path(diag_dir, "unrecognized_activity3_labels.csv")
    )
    
    print(bad_activity_labels)
    stop(
      "Unrecognized activity3 labels found. ",
      "See Diagnostics/unrecognized_activity3_labels.csv. ",
      "Add these labels to state_from_activity3() before generating the figure."
    )
  }
  
  mapped_state_counts <- dt[, .N, by = state][order(state)]
  write_diag_table(
    mapped_state_counts,
    file.path(diag_dir, "mapped_state_counts.csv")
  )
  print(mapped_state_counts)
  
  log_msg(
    "Mapped state counts: ",
    paste(mapped_state_counts$state, mapped_state_counts$N, sep = "=", collapse = "; ")
  )
  
  dt[, day_num_int := parse_day_num(day_num)]
  
  bad_daynum <- dt[is.na(day_num_int), .N]
  if (bad_daynum > 0) {
    log_msg("Rows with invalid day_num removed: ", bad_daynum)
  }
  dt <- dt[!is.na(day_num_int)]
  
  # Do not show baseline where raw HR is missing.
  dt[!is.finite(raw_hr), baseline_hr := NA_real_]
  
  # Time within study-day slot.
  dt[, sec_of_day := hour(dt_time) * 3600L + minute(dt_time) * 60L + second(dt_time)]
  dt[, sec_of_day := as.integer(sec_of_day)]
  
  # Clamp edge cases to the last bin of the day.
  dt[sec_of_day >= 24 * 3600, sec_of_day := 24 * 3600 - RES_SECONDS]
  dt[sec_of_day < 0, sec_of_day := 0L]
  
  dt[, x_hr := 24 * (day_num_int - 1L) + sec_of_day / 3600]
  
  subj_summary <- dt[, .(
    n_rows = .N,
    n_study_days_present = uniqueN(day_num_int),
    min_day_num = min(day_num_int, na.rm = TRUE),
    max_day_num = max(day_num_int, na.rm = TRUE),
    time_min = min(dt_time, na.rm = TRUE),
    time_max = max(dt_time, na.rm = TRUE),
    n_raw_missing = sum(!is.finite(raw_hr)),
    n_bl_missing = sum(!is.finite(baseline_hr))
  ), by = p_id][order(p_id)]
  
  write_diag_table(
    subj_summary,
    file.path(diag_dir, "subject_summary_before_fill.csv")
  )
  
  dup_check <- dt[, .N, by = .(p_id, day_num_int, sec_of_day)][N > 1]
  write_diag_table(
    dup_check,
    file.path(diag_dir, "duplicate_pid_time_rows.csv")
  )
  log_msg("Duplicate (p_id, day_num_int, sec_of_day) rows: ", nrow(dup_check))
  
  if (nrow(dup_check) > 0) {
    log_msg("Collapsing duplicate (p_id, day_num_int, sec_of_day) rows.")
    
    dt <- dt[, .(
      dt_time = suppressWarnings(min(dt_time, na.rm = TRUE)),
      raw_hr = if (all(!is.finite(raw_hr))) NA_real_ else mean(raw_hr, na.rm = TRUE),
      baseline_hr = if (all(!is.finite(baseline_hr))) NA_real_ else mean(baseline_hr, na.rm = TRUE),
      state = factor(mode_chr(state), levels = STATE_LEVELS),
      day_num = mode_chr(day_num),
      x_hr = mean(x_hr, na.rm = TRUE)
    ), by = .(p_id, day_num_int, sec_of_day)]
    
    dt[!is.finite(raw_hr), baseline_hr := NA_real_]
  }
  
  setorder(dt, p_id, day_num_int, sec_of_day)
  
  # Diagnostics for study-day coverage in slots 1..7.
  subj_summary_slots <- dt[, .(
    n_rows_in_slots = .N,
    days_present = paste(sort(unique(day_num_int)), collapse = ","),
    first_time = min(dt_time, na.rm = TRUE),
    last_time = max(dt_time, na.rm = TRUE)
  ), by = p_id][order(p_id)]
  
  write_diag_table(
    subj_summary_slots,
    file.path(diag_dir, "subject_summary_study_day_slots.csv")
  )
  
  # ----------------------------
  # Baseline diagnostics by STUDY DAY, primary
  # ----------------------------
  baseline_study_day_level <- dt[, .(
    n_rows_observed = .N,
    any_raw_hr = any(is.finite(raw_hr)),
    any_baseline = any(is.finite(baseline_hr)),
    n_nonmiss_baseline_rows = sum(is.finite(baseline_hr)),
    n_unique_baseline_values = uniqueN(baseline_hr[is.finite(baseline_hr)]),
    baseline_values_seen = paste(
      sort(unique(baseline_hr[is.finite(baseline_hr)])),
      collapse = ";"
    )
  ), by = .(p_id, day_num, day_num_int)][order(p_id, day_num_int)]
  
  baseline_study_day_level[, baseline_study_day_missing := any_raw_hr & !any_baseline]
  baseline_study_day_level[, baseline_study_day_inconsistent := n_unique_baseline_values > 1]
  
  write_diag_table(
    baseline_study_day_level,
    file.path(diag_dir, "baseline_subject_study_day_level.csv")
  )
  
  baseline_study_subject_summary <- baseline_study_day_level[, .(
    n_study_days_present = .N,
    n_study_days_with_any_raw_hr = sum(any_raw_hr),
    n_study_days_with_any_baseline = sum(any_baseline),
    n_study_days_missing_baseline_given_raw = sum(baseline_study_day_missing),
    n_study_days_inconsistent_baseline = sum(baseline_study_day_inconsistent),
    pct_study_days_missing_baseline_given_raw =
      round(100 * sum(baseline_study_day_missing) / pmax(sum(any_raw_hr), 1), 2)
  ), by = p_id][order(p_id)]
  
  write_diag_table(
    baseline_study_subject_summary,
    file.path(diag_dir, "baseline_subject_study_day_summary.csv")
  )
  
  baseline_study_overall_summary <- baseline_study_day_level[, .(
    n_subject_study_days = .N,
    n_subject_study_days_with_any_raw_hr = sum(any_raw_hr),
    n_subject_study_days_with_any_baseline = sum(any_baseline),
    n_subject_study_days_missing_baseline_given_raw = sum(baseline_study_day_missing),
    n_subject_study_days_inconsistent_baseline = sum(baseline_study_day_inconsistent),
    pct_subject_study_days_missing_baseline_given_raw =
      round(100 * sum(baseline_study_day_missing) / pmax(sum(any_raw_hr), 1), 2)
  )]
  
  write_diag_table(
    baseline_study_overall_summary,
    file.path(diag_dir, "baseline_study_day_overall_summary.csv")
  )
  
  # ----------------------------
  # Baseline diagnostics by CALENDAR DATE, secondary
  # ----------------------------
  dt[, date_day := as.Date(dt_time, tz = LOCAL_TZ)]
  
  baseline_calendar_day_level <- dt[, .(
    n_rows_observed = .N,
    any_raw_hr = any(is.finite(raw_hr)),
    any_baseline = any(is.finite(baseline_hr)),
    n_nonmiss_baseline_rows = sum(is.finite(baseline_hr)),
    n_unique_baseline_values = uniqueN(baseline_hr[is.finite(baseline_hr)]),
    baseline_values_seen = paste(
      sort(unique(baseline_hr[is.finite(baseline_hr)])),
      collapse = ";"
    )
  ), by = .(p_id, date_day)][order(p_id, date_day)]
  
  baseline_calendar_day_level[, baseline_calendar_day_missing := any_raw_hr & !any_baseline]
  baseline_calendar_day_level[, baseline_calendar_day_inconsistent := n_unique_baseline_values > 1]
  
  write_diag_table(
    baseline_calendar_day_level,
    file.path(diag_dir, "baseline_subject_calendar_day_level.csv")
  )
  
  baseline_calendar_subject_summary <- baseline_calendar_day_level[, .(
    n_calendar_days_present = .N,
    n_calendar_days_with_any_raw_hr = sum(any_raw_hr),
    n_calendar_days_with_any_baseline = sum(any_baseline),
    n_calendar_days_missing_baseline_given_raw = sum(baseline_calendar_day_missing),
    n_calendar_days_inconsistent_baseline = sum(baseline_calendar_day_inconsistent),
    pct_calendar_days_missing_baseline_given_raw =
      round(100 * sum(baseline_calendar_day_missing) / pmax(sum(any_raw_hr), 1), 2)
  ), by = p_id][order(p_id)]
  
  write_diag_table(
    baseline_calendar_subject_summary,
    file.path(diag_dir, "baseline_subject_calendar_day_summary.csv")
  )
  
  baseline_calendar_overall_summary <- baseline_calendar_day_level[, .(
    n_subject_calendar_days = .N,
    n_subject_calendar_days_with_any_raw_hr = sum(any_raw_hr),
    n_subject_calendar_days_with_any_baseline = sum(any_baseline),
    n_subject_calendar_days_missing_baseline_given_raw = sum(baseline_calendar_day_missing),
    n_subject_calendar_days_inconsistent_baseline = sum(baseline_calendar_day_inconsistent),
    pct_subject_calendar_days_missing_baseline_given_raw =
      round(100 * sum(baseline_calendar_day_missing) / pmax(sum(any_raw_hr), 1), 2)
  )]
  
  write_diag_table(
    baseline_calendar_overall_summary,
    file.path(diag_dir, "baseline_calendar_day_overall_summary.csv")
  )
  
  # ----------------------------
  # Fill missing grid within each study-day slot 1..7
  # ----------------------------
  dt_filled <- rbindlist(
    lapply(split(dt, dt$p_id), fill_missing_grid_one_subject, step_secs = RES_SECONDS),
    use.names = TRUE,
    fill = TRUE
  )
  
  setorder(dt_filled, p_id, day_num_int, sec_of_day)
  dt_filled[, x_hr := 24 * (day_num_int - 1L) + sec_of_day / 3600]
  
  fill_diag <- dt_filled[, .(
    n_total_points = .N,
    n_missing_raw = sum(!is.finite(raw_hr)),
    pct_missing_raw = round(100 * mean(!is.finite(raw_hr)), 2)
  ), by = p_id][order(p_id)]
  
  write_diag_table(
    fill_diag,
    file.path(diag_dir, "missingness_after_fill.csv")
  )
  
  global_y <- NULL
  if (use_global_y) {
    y_vals <- c(dt_filled$raw_hr, dt_filled$baseline_hr)
    y_vals <- y_vals[is.finite(y_vals)]
    
    if (length(y_vals) > 0) {
      global_y <- pad_range(range(y_vals))
      log_msg("Global y-axis: [", round(global_y[1], 2), ", ", round(global_y[2], 2), "]")
    }
  }
  
  subj_tbl <- data.table(p_id = unique(dt_filled$p_id))
  subj_tbl[, pid_num := pid_sort_key(p_id)]
  setorder(subj_tbl, pid_num, p_id)
  
  subjects <- subj_tbl$p_id
  log_msg("Unique subjects: ", length(subjects))
  
  subject_chunks <- split(subjects, ceiling(seq_along(subjects) / subjects_per_page))
  n_pages <- length(subject_chunks)
  log_msg("Pages to write: ", n_pages)
  
  open_pdf(out_pdf, width = 11, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (pg in seq_along(subject_chunks)) {
    ss <- subject_chunks[[pg]]
    log_msg("Writing page ", pg, " / ", n_pages, " | Subjects: ", paste(ss, collapse = ", "))
    
    header <- make_page_header(pg, n_pages)
    
    plots <- lapply(ss, function(s) {
      make_subject_plot(
        dsub = copy(dt_filled[p_id == s]),
        subject_label = as.character(s),
        global_y = global_y
      )
    })
    
    body <- wrap_plots(plots, ncol = 1)
    print(header / body + plot_layout(heights = c(2.0, 6.5)))
  }
  
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "session_info.txt"))
  log_msg("Session information written to: ", file.path(out_dir, "session_info.txt"))
  log_msg("DONE. PDF written to: ", out_pdf)
}

main()