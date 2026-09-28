###########################################################
# Arctic thermal response curves only
# Cold-side density-adaptive envelope + warm-side fixed 0.1 °C q99 edge
#
# Data source:
# RESULT/Arctic/<Group type>/<Group>/current.img
# EnvironmentData/Arctic/current/bio1.asc
#
# Important differences from the earlier Arctic script:
# 1. Only Arctic groups are calculated and plotted.
# 2. No temperature-tail or suitability outliers are removed.
# 3. Sampling is NOT stratified by suitability class. The default is one
#    random cell per spatial block, which preserves the local suitability
#    density much better than block × suitability-class sampling.
# 4. The cold-side reference curve keeps the density-adaptive local-envelope
#    method: dominant-cluster q99 when suitability is concentrated near zero,
#    all-point q99 when the local distribution is dispersed, and a continuous
#    blend between the two.
# 5. After the cold-side threshold, the warm-side upper edge is calculated
#    independently in fixed 0.1 °C temperature windows. In every window, the
#    99th percentile of ALL suitability values is used; no density-cluster
#    restriction is applied on the warm side.
# 6. The fixed-window warm q99 points are fitted with an anchored monotone
#    Richards curve, with an anchored derivative-integrated spline fallback.
#    This produces one continuous, smooth, non-decreasing upper-edge curve.
# 7. The far cold-side baseline is retained, and the threshold neighbourhood is
#    joined to the warm q99 curve with the existing monotone bridge.
# 8. Threshold extraction is recalculated from the FINAL plotted curve and is
#    strictly the maximum positive normalized curvature on the cold accelerating
#    branch before the maximum-slope point.
# 9. Plot styling parameters are kept identical to the Antarctic panels:
#    grey points, black fitted curve, red long-dashed threshold line, identical
#    point sizes, line widths, axis limits, theme, label size, and output size.
############################################################

# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(terra)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_result_dir <- "E:/1 MY PROJECTS/Climate Change/RESULT"
base_env_dir    <- "E:/1 MY PROJECTS/Climate Change/EnvironmentData"
base_work_dir   <- "E:/1 MY PROJECTS/Climate Change/分析过程"

out_dir <- file.path(
  base_work_dir,
  "_response_curve_figures_Arctic_cold_density_warm_q99"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Groups
# ------------------------------------------------------------

growth_groups <- c("Crustose", "Foliose", "Fruticose")
colour_groups <- c("Pale", "Bright", "Dark")

# ------------------------------------------------------------
# 3. Sampling settings
# ------------------------------------------------------------

# "block_only": one random raster cell per spatial block; recommended because
# it reduces spatial redundancy without distorting the suitability distribution.
# "all_cells": retain every valid raster cell; potentially memory intensive.
arctic_sampling_mode <- "block_only"

block_size_arc <- 5
set.seed(123)

# ------------------------------------------------------------
# 4. Density-adaptive local envelope settings
# ------------------------------------------------------------

# Envelope calculation positions are spaced every 0.1 °C.
temp_step <- 0.10

# Initial moving window is centre ± 0.25 °C, i.e. total width = 0.50 °C.
window_half_width_initial <- 0.25

# If fewer than local_min_n points occur in a window, enlarge the half-width
# incrementally until enough points are available or the maximum is reached.
window_half_width_step <- 0.10
window_half_width_max  <- 1.00
local_min_n            <- 40

# Local upper-envelope quantile.
local_quantile <- 0.99

# The shortest interval containing 70% of local suitability values is used as
# a robust approximation of the dominant high-density cluster.
cluster_mass <- 0.70

# Cluster-width rules on the 0–1 suitability scale:
# width <= 0.20: use the dominant-cluster q99 fully;
# width >= 0.45: use the all-point q99 fully;
# intermediate widths: blend the two estimates continuously.
cluster_width_full <- 0.20
cluster_width_none <- 0.45

# First-stage smoothing of the raw local envelope points.
# This stage follows the local density-adaptive envelope before the monotonic
# constraint is imposed.
envelope_spar <- 0.84

# Cold-side preservation and warm-side smoothing.
# The far-left baseline is copied without modification. Only a narrow interval
# surrounding the threshold is replaced by a monotone C1 Hermite bridge. The
# bridge removes the small hook/corner marked by the arrows while retaining the
# original near-zero cold-side structure outside that narrow interval.
left_preserve_after_threshold     <- 1.00  # °C retained after the threshold
left_preserve_rise_fraction       <- 0.08  # retain at least the first 8% rise
left_preserve_max_after_threshold <- 2.00  # avoid carrying a plateau too far right

# The local C1 bridge begins at the last near-baseline point before the
# threshold, subject to minimum and maximum widths. This adaptive rule captures
# the small hook marked by the arrows without altering the colder flat baseline.
threshold_bridge_baseline_fraction <- 0.003
threshold_bridge_min_left_width     <- 0.70  # °C
threshold_bridge_max_left_width     <- 2.20  # °C

# Number of grid intervals used to estimate the endpoint slopes of the bridge.
threshold_bridge_slope_span <- 4L

# Growth-form panels retain a slightly stronger residual corner at the cold-side
# departure, especially for Crustose. For Growth form only, start the bridge
# farther into the exact near-zero baseline and force a horizontal cold-end
# tangent. Color-type curve geometry keeps the previous settings unchanged.
growth_bridge_baseline_fraction <- 0.0015
growth_bridge_min_left_width     <- 1.40  # °C
growth_bridge_max_left_width     <- 3.00  # °C
growth_bridge_force_zero_left_slope <- TRUE

# Warm-side upper-edge calculation:
# centres are separated by 0.1 °C and each fixed window has total width 0.1 °C.
# The target is the 99th percentile of ALL suitability values in that window,
# i.e. the lower boundary of the locally highest 1% suitability values.
warm_bin_width        <- 0.10
warm_bin_half_width   <- warm_bin_width / 2
warm_min_bin_n        <- 8L
warm_quantile         <- 0.99

# The warm q99 fit starts at the first 0.1 °C position after the preliminary
# cold-side threshold. This makes the post-threshold curve directly controlled
# by the fixed-window all-point q99 values rather than by the cold reference.
warm_start_after_threshold <- warm_bin_width

# Diagnostic statistic only: mean of the highest 1% values in each fixed bin.
# The fitted upper edge itself uses q99_all, not this mean.
warm_top_fraction     <- 0.01

# Warm q99 curve fitting. The anchored Richards curve is monotone and flexible:
# shape = 1 gives an early concave-down rise, while shape > 1 permits a smoother
# accelerating transition before saturation.
right_curve_method          <- "anchored_richards_q99"
right_min_fit_span          <- 4.00
right_anchor_weight         <- 800
warm_underfit_penalty       <- 1.75
warm_richards_rate_bounds   <- c(0.15, 20.0)
warm_richards_shape_bounds  <- c(1.00, 6.00)

# Fallback settings if the Richards optimization fails.
warm_spline_spar            <- 0.86
warm_derivative_spar        <- 0.90

# Robust warm-end anchor based only on the fixed-window q99 points.
warm_upper_anchor_fraction  <- 0.10
warm_upper_anchor_quantile  <- 0.75

# Number of points in the final fitted curve.
grid_n <- 800

# ------------------------------------------------------------
# 5. Threshold settings
# Same settings as the revised Antarctic cold-side method
# ------------------------------------------------------------

threshold_edge_trim         <- 0.02
threshold_rise_trim         <- c(0.005, 0.30)
threshold_min_slope_fraction <- 0.03

# ------------------------------------------------------------
# 6. Plot settings
# Identical to the Antarctic figure settings
# ------------------------------------------------------------

max_scatter_n <- 5000

point_col     <- "grey55"
curve_col     <- "black"
threshold_col <- "#FF4D4D"

point_size       <- 0.85
point_alpha      <- 0.32
curve_linewidth  <- 1.10
threshold_width  <- 0.85
threshold_type   <- "longdash"
threshold_alpha  <- 0.95

figure_width_mm  <- 180
figure_height_mm <- 85
figure_dpi       <- 300

# ------------------------------------------------------------
# 7. Path helpers
# ------------------------------------------------------------

make_suitability_path <- function(group_type_name, group_name) {
  file.path(
    base_result_dir,
    "Arctic",
    group_type_name,
    group_name,
    "current.img"
  )
}

make_bio1_path <- function() {
  file.path(
    base_env_dir,
    "Arctic",
    "current",
    "bio1.asc"
  )
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File not found:\n", path, call. = FALSE)
  }
}

read_project_to_template <- function(r_path, template_r, method = "bilinear") {
  check_file_exists(r_path)
  r <- terra::rast(r_path)
  terra::project(r, template_r, method = method)
}

# ------------------------------------------------------------
# 8. Build Arctic sample
# No suitability-class stratification and no outlier deletion
# ------------------------------------------------------------

build_arctic_sample <- function(group_type_name, group_name) {
  su_path   <- make_suitability_path(group_type_name, group_name)
  bio1_path <- make_bio1_path()
  
  check_file_exists(su_path)
  check_file_exists(bio1_path)
  
  su_current <- terra::rast(su_path)
  bio1_r <- read_project_to_template(
    r_path = bio1_path,
    template_r = su_current,
    method = "bilinear"
  )
  
  names(su_current) <- "suitability_raw"
  names(bio1_r) <- "bio1"
  
  su_values <- terra::values(su_current, mat = FALSE)
  cell_id <- which(!is.na(su_values))
  
  if (length(cell_id) == 0) {
    stop("No valid suitability cells found: ", su_path, call. = FALSE)
  }
  
  xy <- terra::xyFromCell(su_current, cell_id)
  
  su_df <- tibble(
    cell = cell_id,
    x = xy[, 1],
    y = xy[, 2],
    suitability_raw = as.numeric(su_values[cell_id])
  ) %>%
    filter(
      !is.na(suitability_raw),
      is.finite(suitability_raw),
      suitability_raw >= 0,
      suitability_raw <= 1000
    ) %>%
    mutate(
      suit01 = suitability_raw / 1000,
      suit01_clip = pmin(pmax(suit01, 0.001), 0.999)
    )
  
  rc <- terra::rowColFromCell(su_current, su_df$cell)
  
  su_df <- su_df %>%
    mutate(
      row = rc[, 1],
      col = rc[, 2],
      block_row = floor((row - 1) / block_size_arc),
      block_col = floor((col - 1) / block_size_arc),
      block_id = paste(block_row, block_col, sep = "_")
    )
  
  if (identical(arctic_sampling_mode, "block_only")) {
    sample_df <- su_df %>%
      group_by(block_id) %>%
      slice_sample(n = 1) %>%
      ungroup()
  } else if (identical(arctic_sampling_mode, "all_cells")) {
    sample_df <- su_df
  } else {
    stop(
      "Unknown arctic_sampling_mode: ", arctic_sampling_mode,
      ". Use 'block_only' or 'all_cells'.",
      call. = FALSE
    )
  }
  
  pts <- terra::vect(
    sample_df,
    geom = c("x", "y"),
    crs = terra::crs(su_current)
  )
  
  bio1_vals <- terra::extract(bio1_r, pts)[, 2]
  
  sample_out <- sample_df %>%
    mutate(bio1 = as.numeric(bio1_vals)) %>%
    filter(
      !is.na(bio1),
      is.finite(bio1),
      bio1 > -9000,
      bio1 < 10000
    )
  
  # WorldClim bio1 may be stored as °C × 10.
  if (median(abs(sample_out$bio1), na.rm = TRUE) > 80) {
    message(
      "Converting bio1 from °C × 10 to °C for Arctic | ",
      group_type_name, " | ", group_name
    )
    sample_out <- sample_out %>% mutate(bio1 = bio1 / 10)
  }
  
  message(
    "\nPrepared Arctic sample: ", group_type_name, " | ", group_name,
    "\n  sampling mode = ", arctic_sampling_mode,
    "\n  n = ", nrow(sample_out),
    "\n  suitability range = ",
    paste(round(range(sample_out$suit01, na.rm = TRUE), 3), collapse = " to "),
    "\n  bio1 range = ",
    paste(round(range(sample_out$bio1, na.rm = TRUE), 3), collapse = " to ")
  )
  
  sample_out
}

# ------------------------------------------------------------
# 9. Local-density helpers
# ------------------------------------------------------------

# Return the shortest interval containing the requested fraction of values.
# This is a robust one-dimensional high-density interval approximation.
shortest_mass_interval <- function(y, mass = cluster_mass) {
  y <- sort(y[is.finite(y)])
  n <- length(y)
  
  if (n < 2) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      width = NA_real_,
      n_interval = 0L
    ))
  }
  
  m <- ceiling(mass * n)
  m <- min(max(m, 2L), n)
  
  start_index <- seq_len(n - m + 1L)
  end_index <- start_index + m - 1L
  widths <- y[end_index] - y[start_index]
  
  best <- which.min(widths)
  
  list(
    lower = y[start_index[best]],
    upper = y[end_index[best]],
    width = widths[best],
    n_interval = m
  )
}

# Convert dominant-cluster width into a continuous cluster weight.
# 1 = completely use cluster q99; 0 = completely use all-point q99.
cluster_weight_from_width <- function(width) {
  if (!is.finite(width)) return(0)
  
  if (width <= cluster_width_full) return(1)
  if (width >= cluster_width_none) return(0)
  
  1 - (width - cluster_width_full) /
    (cluster_width_none - cluster_width_full)
}

# Extract the local values around one temperature centre, expanding the window
# only when needed to reach the minimum sample size.
get_local_window <- function(dat, centre) {
  half_width <- window_half_width_initial
  
  repeat {
    idx <- abs(dat$bio1 - centre) <= half_width
    n_here <- sum(idx)
    
    if (n_here >= local_min_n || half_width >= window_half_width_max) {
      break
    }
    
    half_width <- min(
      half_width + window_half_width_step,
      window_half_width_max
    )
  }
  
  list(
    data = dat[idx, , drop = FALSE],
    half_width = half_width
  )
}

calculate_local_envelope_points <- function(sample_out) {
  dat <- sample_out %>%
    transmute(
      bio1 = as.numeric(bio1),
      suit01 = pmin(pmax(as.numeric(suit01), 0), 1)
    ) %>%
    filter(
      is.finite(bio1),
      is.finite(suit01)
    ) %>%
    arrange(bio1)
  
  if (nrow(dat) < local_min_n) {
    stop(
      "Too few Arctic sample points for local-envelope calculation: n = ",
      nrow(dat), ".",
      call. = FALSE
    )
  }
  
  temp_min <- ceiling(min(dat$bio1, na.rm = TRUE) / temp_step) * temp_step
  temp_max <- floor(max(dat$bio1, na.rm = TRUE) / temp_step) * temp_step
  
  if (!is.finite(temp_min) || !is.finite(temp_max) || temp_min >= temp_max) {
    stop("Invalid Arctic temperature range.", call. = FALSE)
  }
  
  centres <- seq(temp_min, temp_max, by = temp_step)
  result <- vector("list", length(centres))
  
  for (i in seq_along(centres)) {
    centre <- centres[i]
    local_obj <- get_local_window(dat, centre)
    local_dat <- local_obj$data
    
    if (nrow(local_dat) < max(10L, ceiling(local_min_n / 2))) {
      next
    }
    
    y <- local_dat$suit01
    interval_obj <- shortest_mass_interval(y, mass = cluster_mass)
    
    cluster_values <- y[
      y >= interval_obj$lower &
        y <= interval_obj$upper
    ]
    
    q_all <- as.numeric(
      quantile(y, probs = local_quantile, na.rm = TRUE, type = 8)
    )
    
    q_cluster <- if (length(cluster_values) >= 5) {
      as.numeric(
        quantile(
          cluster_values,
          probs = local_quantile,
          na.rm = TRUE,
          type = 8
        )
      )
    } else {
      q_all
    }
    
    cluster_weight <- cluster_weight_from_width(interval_obj$width)
    
    # Smooth transition between the dominant-cluster envelope and the complete
    # local upper envelope. This avoids an artificial jump when the local
    # distribution changes gradually from concentrated to dispersed.
    envelope_value <-
      cluster_weight * q_cluster +
      (1 - cluster_weight) * q_all
    
    result[[i]] <- tibble(
      bio1 = centre,
      envelope_raw = pmin(pmax(envelope_value, 0), 1),
      q99_cluster = pmin(pmax(q_cluster, 0), 1),
      q99_all = pmin(pmax(q_all, 0), 1),
      cluster_weight = cluster_weight,
      clustered = cluster_weight >= 0.5,
      cluster_lower = interval_obj$lower,
      cluster_upper = interval_obj$upper,
      cluster_width = interval_obj$width,
      n_window = nrow(local_dat),
      n_cluster = length(cluster_values),
      window_half_width = local_obj$half_width
    )
  }
  
  envelope_points <- bind_rows(result) %>%
    filter(
      is.finite(bio1),
      is.finite(envelope_raw)
    ) %>%
    arrange(bio1)
  
  if (nrow(envelope_points) < 8) {
    stop(
      "Too few usable local-envelope points: n = ",
      nrow(envelope_points), ".",
      call. = FALSE
    )
  }
  
  envelope_points
}

# Fixed 0.1 °C warm-side windows. Unlike the cold-side density-adaptive
# calculation, this function always uses ALL suitability values in each bin.
# The q99 value is the lower boundary of the locally highest 1% values.
calculate_fixed_warm_q99_points <- function(
    sample_out,
    threshold_x = -Inf
) {
  dat <- sample_out %>%
    transmute(
      bio1 = as.numeric(bio1),
      suit01 = pmin(pmax(as.numeric(suit01), 0), 1)
    ) %>%
    filter(
      is.finite(bio1),
      is.finite(suit01)
    ) %>%
    arrange(bio1)
  
  if (nrow(dat) < warm_min_bin_n) {
    stop(
      "Too few Arctic sample points for warm q99 calculation: n = ",
      nrow(dat), ".",
      call. = FALSE
    )
  }
  
  centre_min <- ceiling(
    (min(dat$bio1, na.rm = TRUE) + warm_bin_half_width) /
      warm_bin_width
  ) * warm_bin_width
  
  centre_max <- floor(
    (max(dat$bio1, na.rm = TRUE) - warm_bin_half_width) /
      warm_bin_width
  ) * warm_bin_width
  
  centres <- seq(centre_min, centre_max, by = warm_bin_width)
  result <- vector("list", length(centres))
  
  for (i in seq_along(centres)) {
    centre <- centres[i]
    
    # Half-open bins prevent a point exactly on a boundary from being counted
    # in two adjacent 0.1 °C windows.
    idx <- dat$bio1 >= (centre - warm_bin_half_width) &
      dat$bio1 < (centre + warm_bin_half_width)
    
    local_y <- dat$suit01[idx]
    local_y <- local_y[is.finite(local_y)]
    n_here <- length(local_y)
    
    if (n_here < warm_min_bin_n) {
      next
    }
    
    q99_all <- as.numeric(
      quantile(
        local_y,
        probs = warm_quantile,
        na.rm = TRUE,
        type = 8
      )
    )
    
    n_top <- max(1L, ceiling(warm_top_fraction * n_here))
    top_values <- sort(local_y, decreasing = TRUE)[seq_len(n_top)]
    
    result[[i]] <- tibble(
      bio1 = centre,
      q99_all = pmin(pmax(q99_all, 0), 1),
      top1_mean = pmin(pmax(mean(top_values, na.rm = TRUE), 0), 1),
      top1_min = pmin(pmax(min(top_values, na.rm = TRUE), 0), 1),
      top1_max = pmin(pmax(max(top_values, na.rm = TRUE), 0), 1),
      n_bin = n_here,
      bin_lower = centre - warm_bin_half_width,
      bin_upper = centre + warm_bin_half_width
    )
  }
  
  warm_points <- bind_rows(result) %>%
    filter(
      is.finite(bio1),
      is.finite(q99_all),
      bio1 >= threshold_x - warm_bin_width
    ) %>%
    arrange(bio1)
  
  if (nrow(warm_points) < 6) {
    stop(
      "Too few fixed 0.1 °C warm q99 points after the threshold: n = ",
      nrow(warm_points), ". Consider reducing warm_min_bin_n.",
      call. = FALSE
    )
  }
  
  warm_points
}

# ------------------------------------------------------------
# 10. Preserve the original cold side; fit warm fixed-window q99 edge
# ------------------------------------------------------------

# Anchored Richards transformation. It is exactly y0 at x0 and exactly
# upper at x1. The rate controls how quickly the curve approaches the upper
# edge; shape controls whether the rise begins immediately or accelerates first.
predict_anchored_richards <- function(
    x,
    upper,
    rate,
    shape,
    x0,
    x1,
    y0
) {
  rate <- max(as.numeric(rate), 1e-6)
  shape <- max(as.numeric(shape), 1e-6)
  
  t <- (x - x0) / max(x1 - x0, 1e-8)
  t <- pmin(pmax(t, 0), 1)
  
  numerator <- (1 - exp(-rate * t))^shape
  denominator <- (1 - exp(-rate))^shape
  
  if (!is.finite(denominator) || denominator <= 1e-12) {
    g <- t
  } else {
    g <- numerator / denominator
  }
  
  g <- pmin(pmax(g, 0), 1)
  y0 + (upper - y0) * g
}

# Smoothstep blend: weight and its first derivative are both zero at the start,
# so the preserved reference curve joins the warm-side curve continuously.
smoothstep_weight <- function(x, start, end) {
  if (!is.finite(end) || end <= start) {
    return(as.numeric(x > start))
  }
  
  t <- (x - start) / (end - start)
  t <- pmin(pmax(t, 0), 1)
  t^2 * (3 - 2 * t)
}

# Estimate a local non-negative slope at one curve index.
estimate_local_nonnegative_slope <- function(
    x,
    y,
    index,
    span = threshold_bridge_slope_span,
    direction = c("both", "backward", "forward")
) {
  direction <- match.arg(direction)
  n <- length(x)
  index <- max(1L, min(as.integer(index), n))
  span <- max(1L, as.integer(span))
  
  if (direction == "backward") {
    i0 <- max(1L, index - span)
    i1 <- index
  } else if (direction == "forward") {
    i0 <- index
    i1 <- min(n, index + span)
  } else {
    i0 <- max(1L, index - span)
    i1 <- min(n, index + span)
  }
  
  if (i1 <= i0 || !is.finite(x[i1] - x[i0]) || x[i1] <= x[i0]) {
    return(0)
  }
  
  slope <- (y[i1] - y[i0]) / (x[i1] - x[i0])
  if (!is.finite(slope)) slope <- 0
  max(slope, 0)
}

# Monotone cubic Hermite bridge. Endpoint slopes are Fritsch-Carlson limited,
# which prevents overshoot and guarantees a smooth non-decreasing connection.
predict_monotone_hermite_bridge <- function(
    x,
    x0,
    x1,
    y0,
    y1,
    m0,
    m1
) {
  if (!is.finite(x1 - x0) || x1 <= x0 || !is.finite(y1 - y0)) {
    return(rep(y0, length(x)))
  }
  
  interval <- x1 - x0
  delta <- (y1 - y0) / interval
  
  if (!is.finite(delta) || delta <= 0) {
    t <- pmin(pmax((x - x0) / interval, 0), 1)
    return(y0 + (y1 - y0) * t)
  }
  
  m0 <- max(as.numeric(m0), 0)
  m1 <- max(as.numeric(m1), 0)
  
  # Fritsch-Carlson monotonicity limiter.
  alpha <- m0 / delta
  beta  <- m1 / delta
  norm2 <- alpha^2 + beta^2
  
  if (is.finite(norm2) && norm2 > 9) {
    tau <- 3 / sqrt(norm2)
    m0 <- tau * alpha * delta
    m1 <- tau * beta  * delta
  }
  
  m0 <- min(m0, 3 * delta)
  m1 <- min(m1, 3 * delta)
  
  t <- pmin(pmax((x - x0) / interval, 0), 1)
  t2 <- t^2
  t3 <- t^3
  
  h00 <-  2 * t3 - 3 * t2 + 1
  h10 <-      t3 - 2 * t2 + t
  h01 <- -2 * t3 + 3 * t2
  h11 <-      t3 -     t2
  
  y <- h00 * y0 + h10 * interval * m0 +
    h01 * y1 + h11 * interval * m1
  
  pmin(pmax(y, min(y0, y1)), max(y0, y1))
}

# Quintic Hermite bridge used for Growth form. With a zero cold-end slope and
# zero endpoint second derivatives, it joins the exact flat baseline without
# the small residual corner that remained in the Crustose panel. Slopes are
# limited relative to the secant slope and a final monotone projection is used
# only inside the bridge interval.
predict_c2_quintic_bridge <- function(
    x,
    x0,
    x1,
    y0,
    y1,
    m0,
    m1
) {
  if (!is.finite(x1 - x0) || x1 <= x0 || !is.finite(y1 - y0)) {
    return(rep(y0, length(x)))
  }
  
  h <- x1 - x0
  delta <- (y1 - y0) / h
  
  if (!is.finite(delta) || delta <= 0) {
    t <- pmin(pmax((x - x0) / h, 0), 1)
    return(y0 + (y1 - y0) * (6 * t^5 - 15 * t^4 + 10 * t^3))
  }
  
  m0 <- pmin(pmax(as.numeric(m0), 0), 2.5 * delta)
  m1 <- pmin(pmax(as.numeric(m1), 0), 2.5 * delta)
  
  t <- pmin(pmax((x - x0) / h, 0), 1)
  t2 <- t^2
  t3 <- t^3
  t4 <- t^4
  t5 <- t^5
  
  h00 <- 1 - 10 * t3 + 15 * t4 - 6 * t5
  h10 <- t - 6 * t3 + 8 * t4 - 3 * t5
  h01 <- 10 * t3 - 15 * t4 + 6 * t5
  h11 <- -4 * t3 + 7 * t4 - 3 * t5
  
  y <- h00 * y0 + h10 * h * m0 +
    h01 * y1 + h11 * h * m1
  
  y <- pmin(pmax(y, min(y0, y1)), max(y0, y1))
  y <- cummax(y)
  
  # Preserve exact endpoint values after monotone projection.
  if (length(y) >= 1) y[1] <- y0
  if (length(y) >= 2) y[length(y)] <- y1
  
  y
}

# Fit the warm-side segment to the fixed 0.1 °C all-point q99 values.
# An anchored Richards curve is preferred. The fallback smooths the q99 points,
# extracts a non-negative smooth derivative, and integrates that derivative to
# obtain a continuous monotone curve with exact endpoint anchors.
fit_warm_q99_curve <- function(
    warm_q99_points,
    x_seq,
    splice_x,
    y_join,
    upper_anchor
) {
  x1 <- max(x_seq, na.rm = TRUE)
  
  obs <- warm_q99_points %>%
    filter(bio1 >= splice_x) %>%
    transmute(
      x = as.numeric(bio1),
      y = pmin(pmax(as.numeric(q99_all), y_join), 1),
      w = sqrt(pmax(as.numeric(n_bin), 1))
    ) %>%
    filter(
      is.finite(x),
      is.finite(y),
      is.finite(w),
      w > 0
    ) %>%
    arrange(x)
  
  if (nrow(obs) < 6 || (x1 - splice_x) < 0.5) {
    g <- (x_seq - splice_x) / max(x1 - splice_x, 1e-8)
    g <- pmin(pmax(g, 0), 1)
    
    return(list(
      y = y_join + (upper_anchor - y_join) * g,
      method = "linear anchored fixed-window q99 fallback",
      parameters = c(upper = upper_anchor, rate = NA, shape = NA)
    ))
  }
  
  if (identical(right_curve_method, "anchored_richards_q99") &&
      upper_anchor > y_join + 1e-4) {
    
    upper_start <- min(
      1,
      max(
        upper_anchor,
        as.numeric(quantile(obs$y, 0.85, na.rm = TRUE, type = 8))
      )
    )
    
    upper_lower <- min(
      0.999,
      max(
        y_join + 0.01,
        as.numeric(quantile(obs$y, 0.60, na.rm = TRUE, type = 8))
      )
    )
    
    objective <- function(par) {
      pred <- predict_anchored_richards(
        x = obs$x,
        upper = par[1],
        rate = par[2],
        shape = par[3],
        x0 = splice_x,
        x1 = x1,
        y0 = y_join
      )
      
      residual <- obs$y - pred
      
      # A q99 curve is an upper-edge target. Underprediction is therefore given
      # a moderately larger penalty than overprediction, without forcing the
      # curve through isolated single-bin maxima.
      penalty <- ifelse(
        residual > 0,
        warm_underfit_penalty,
        1
      )
      
      sum(obs$w * penalty * residual^2, na.rm = TRUE) /
        sum(obs$w, na.rm = TRUE)
    }
    
    opt <- tryCatch(
      optim(
        par = c(
          upper_start,
          4.0,
          1.5
        ),
        fn = objective,
        method = "L-BFGS-B",
        lower = c(
          upper_lower,
          warm_richards_rate_bounds[1],
          warm_richards_shape_bounds[1]
        ),
        upper = c(
          1,
          warm_richards_rate_bounds[2],
          warm_richards_shape_bounds[2]
        )
      ),
      error = function(e) NULL
    )
    
    if (!is.null(opt) &&
        is.finite(opt$value) &&
        opt$convergence == 0) {
      
      y_richards <- predict_anchored_richards(
        x = x_seq,
        upper = opt$par[1],
        rate = opt$par[2],
        shape = opt$par[3],
        x0 = splice_x,
        x1 = x1,
        y0 = y_join
      )
      
      return(list(
        y = pmin(pmax(y_richards, y_join), 1),
        method = "anchored Richards fit to fixed 0.1 °C all-point q99",
        parameters = c(
          upper = opt$par[1],
          rate = opt$par[2],
          shape = opt$par[3]
        )
      ))
    }
  }
  
  # Fallback: smooth the fixed-window q99 points, then reconstruct the curve
  # from a smoothed non-negative derivative. This avoids the steps produced by
  # direct isotonic regression while guaranteeing a continuous monotone rise.
  fit_tbl <- bind_rows(
    tibble(x = splice_x, y = y_join, w = right_anchor_weight),
    obs,
    tibble(x = x1, y = upper_anchor, w = right_anchor_weight)
  ) %>%
    group_by(x) %>%
    summarise(
      y = weighted.mean(y, w, na.rm = TRUE),
      w = sum(w, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(x)
  
  sp_level <- smooth.spline(
    x = fit_tbl$x,
    y = fit_tbl$y,
    w = fit_tbl$w,
    spar = warm_spline_spar
  )
  
  warm_idx <- x_seq >= splice_x
  x_warm <- x_seq[warm_idx]
  
  d_raw <- as.numeric(
    predict(sp_level, x = x_warm, deriv = 1)$y
  )
  d_raw[!is.finite(d_raw)] <- 0
  d_raw <- pmax(d_raw, 0)
  
  sp_derivative <- tryCatch(
    smooth.spline(
      x = x_warm,
      y = d_raw,
      spar = warm_derivative_spar
    ),
    error = function(e) NULL
  )
  
  d_smooth <- if (is.null(sp_derivative)) {
    d_raw
  } else {
    as.numeric(predict(sp_derivative, x = x_warm)$y)
  }
  
  d_smooth[!is.finite(d_smooth)] <- 0
  d_smooth <- pmax(d_smooth, 0)
  
  dx <- diff(x_warm)
  cumulative_rise <- c(
    0,
    cumsum(
      0.5 * (d_smooth[-1] + d_smooth[-length(d_smooth)]) * dx
    )
  )
  
  if (!is.finite(max(cumulative_rise, na.rm = TRUE)) ||
      max(cumulative_rise, na.rm = TRUE) <= 1e-12) {
    cumulative_rise <- seq(0, 1, length.out = length(x_warm))
  } else {
    cumulative_rise <- cumulative_rise /
      max(cumulative_rise, na.rm = TRUE)
  }
  
  y_warm <- y_join +
    (upper_anchor - y_join) * cumulative_rise
  
  y_out <- rep(y_join, length(x_seq))
  y_out[warm_idx] <- y_warm
  
  list(
    y = pmin(pmax(y_out, y_join), 1),
    method = paste0(
      "anchored derivative-integrated spline fit to fixed 0.1 °C ",
      "all-point q99"
    ),
    parameters = c(
      upper = upper_anchor,
      rate = NA,
      shape = NA
    )
  )
}

fit_density_adaptive_envelope <- function(
    sample_out,
    group_type_name,
    group_name
) {
  is_growth_form <- identical(group_type_name, "Growth form")
  
  bridge_baseline_fraction_use <- if (is_growth_form) {
    growth_bridge_baseline_fraction
  } else {
    threshold_bridge_baseline_fraction
  }
  
  bridge_min_left_width_use <- if (is_growth_form) {
    growth_bridge_min_left_width
  } else {
    threshold_bridge_min_left_width
  }
  
  bridge_max_left_width_use <- if (is_growth_form) {
    growth_bridge_max_left_width
  } else {
    threshold_bridge_max_left_width
  }
  
  envelope_points <- calculate_local_envelope_points(sample_out)
  
  # Stage 1 is exactly the same as in the previous density-adaptive local 99%
  # envelope code. This creates the reference curve whose cold side is retained.
  fit_weights <- sqrt(pmax(envelope_points$n_window, 1))
  
  smooth_fit <- tryCatch(
    smooth.spline(
      x = envelope_points$bio1,
      y = envelope_points$envelope_raw,
      w = fit_weights,
      spar = envelope_spar
    ),
    error = function(e) {
      stop(
        "Local-envelope smoothing failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  envelope_points$envelope_smooth <- pmin(
    pmax(
      as.numeric(
        predict(smooth_fit, x = envelope_points$bio1)$y
      ),
      0
    ),
    1
  )
  
  iso_fit <- isoreg(
    envelope_points$bio1,
    envelope_points$envelope_smooth
  )
  
  envelope_points$envelope_reference_points <- pmin(
    pmax(as.numeric(iso_fit$yf), 0),
    1
  )
  
  reference_fun <- splinefun(
    x = envelope_points$bio1,
    y = envelope_points$envelope_reference_points,
    method = "monoH.FC"
  )
  
  x_seq <- seq(
    min(envelope_points$bio1, na.rm = TRUE),
    max(envelope_points$bio1, na.rm = TRUE),
    length.out = grid_n
  )
  
  y_reference <- pmin(
    pmax(as.numeric(reference_fun(x_seq)), 0),
    1
  )
  
  reference_curve_df <- tibble(
    bio1 = x_seq,
    suit01 = y_reference
  )
  
  # The threshold is calculated from the original reference curve before any
  # warm-side replacement. Therefore both the threshold and its blue-box cold
  # segment are protected from later smoothing.
  # Preliminary curvature threshold on the unmodified reference curve. It is
  # used only to position the local transition bridge. The reported threshold
  # is recalculated later from the final plotted curve.
  threshold_reference <- find_coldside_curvature_threshold(
    reference_curve_df
  )
  
  threshold_seed_x <- threshold_reference$threshold
  threshold_x <- threshold_seed_x
  
  y_ref_min <- min(y_reference, na.rm = TRUE)
  y_ref_max <- max(y_reference, na.rm = TRUE)
  y_ref_rng <- y_ref_max - y_ref_min
  
  rise_fraction <- if (is.finite(y_ref_rng) && y_ref_rng > 0) {
    (y_reference - y_ref_min) / y_ref_rng
  } else {
    rep(0, length(y_reference))
  }
  
  # The post-threshold upper edge must be controlled by the fixed 0.1 °C
  # all-point q99 values. Therefore the warm fit starts at the first 0.1 °C
  # position after the preliminary cold-side threshold, rather than one or two
  # degrees farther to the right as in the earlier preservation script.
  if (is.finite(threshold_x)) {
    splice_target <- threshold_x + warm_start_after_threshold
  } else {
    rise_candidates <- which(
      rise_fraction >= left_preserve_rise_fraction
    )
    
    splice_target <- if (length(rise_candidates) > 0) {
      x_seq[rise_candidates[1]]
    } else {
      min(x_seq) + warm_start_after_threshold
    }
  }
  
  # Leave enough temperature span for the warm-side fit.
  splice_target <- min(
    splice_target,
    max(x_seq) - right_min_fit_span
  )
  
  splice_index <- max(which(x_seq <= splice_target))
  splice_index <- max(2L, min(splice_index, length(x_seq) - 2L))
  splice_x <- x_seq[splice_index]
  y_join <- y_reference[splice_index]
  
  # Fixed 0.1 °C all-point q99 values are calculated independently of the
  # cold-side density-adaptive envelope. Only bins after the threshold are used
  # to define the warm upper edge.
  warm_q99_points <- calculate_fixed_warm_q99_points(
    sample_out = sample_out,
    threshold_x = threshold_seed_x
  )
  
  warm_q99_fit_points <- warm_q99_points %>%
    filter(bio1 >= splice_x)
  
  if (nrow(warm_q99_fit_points) < 6) {
    stop(
      "Too few warm q99 bins at or after the splice: n = ",
      nrow(warm_q99_fit_points), ".",
      call. = FALSE
    )
  }
  
  n_upper <- max(
    5L,
    floor(
      warm_upper_anchor_fraction *
        nrow(warm_q99_fit_points)
    )
  )
  
  upper_anchor <- as.numeric(
    quantile(
      tail(warm_q99_fit_points$q99_all, n_upper),
      probs = warm_upper_anchor_quantile,
      na.rm = TRUE,
      type = 8
    )
  )
  
  upper_anchor <- pmin(
    pmax(upper_anchor, y_join + 0.02),
    1
  )
  
  warm_fit <- fit_warm_q99_curve(
    warm_q99_points = warm_q99_points,
    x_seq = x_seq,
    splice_x = splice_x,
    y_join = y_join,
    upper_anchor = upper_anchor
  )
  
  # The previous version copied the threshold neighbourhood exactly. That also
  # retained a small angular kink caused by the isotonic reference curve. Here,
  # only the far cold-side baseline is copied exactly. A narrow monotone C1
  # Hermite bridge replaces the threshold neighbourhood and joins the anchored
  # warm-side curve with matched value and slope.
  bridge_start_target <- if (is.finite(threshold_x)) {
    baseline_candidates <- which(
      x_seq <= threshold_x &
        rise_fraction <= bridge_baseline_fraction_use
    )
    
    baseline_start <- if (length(baseline_candidates) > 0) {
      x_seq[max(baseline_candidates)]
    } else {
      threshold_x - bridge_min_left_width_use
    }
    
    # Use the colder of the baseline-derived point and the minimum-width point,
    # but do not let the bridge extend excessively far into the cold tail.
    max(
      min(
        baseline_start,
        threshold_x - bridge_min_left_width_use
      ),
      threshold_x - bridge_max_left_width_use
    )
  } else {
    splice_x - bridge_min_left_width_use
  }
  
  bridge_start_index <- max(which(x_seq <= bridge_start_target))
  bridge_start_index <- max(
    2L,
    min(bridge_start_index, splice_index - 1L)
  )
  
  bridge_start_x <- x_seq[bridge_start_index]
  bridge_start_y <- y_reference[bridge_start_index]
  
  bridge_left_slope_estimated <- estimate_local_nonnegative_slope(
    x = x_seq,
    y = y_reference,
    index = bridge_start_index,
    direction = "backward"
  )
  
  # The small arrow-marked hook in the Growth-form panels is caused by a
  # non-zero inherited tangent at the exact-baseline/bridge junction. A zero
  # cold-end tangent removes that residual corner while the colder baseline is
  # still copied bit-for-bit. Color-type curves retain the previous estimate.
  bridge_left_slope <- if (
    is_growth_form && growth_bridge_force_zero_left_slope
  ) {
    0
  } else {
    bridge_left_slope_estimated
  }
  
  bridge_right_slope <- estimate_local_nonnegative_slope(
    x = x_seq,
    y = warm_fit$y,
    index = splice_index,
    direction = "forward"
  )
  
  bridge_idx <- bridge_start_index:splice_index
  
  if (is_growth_form) {
    bridge_y <- predict_c2_quintic_bridge(
      x = x_seq[bridge_idx],
      x0 = bridge_start_x,
      x1 = splice_x,
      y0 = bridge_start_y,
      y1 = y_join,
      m0 = bridge_left_slope,
      m1 = bridge_right_slope
    )
    bridge_method <- "growth-form zero-slope C2 quintic bridge"
  } else {
    bridge_y <- predict_monotone_hermite_bridge(
      x = x_seq[bridge_idx],
      x0 = bridge_start_x,
      x1 = splice_x,
      y0 = bridge_start_y,
      y1 = y_join,
      m0 = bridge_left_slope,
      m1 = bridge_right_slope
    )
    bridge_method <- "monotone C1 threshold bridge"
  }
  
  # Exact-copy rule now applies to the far cold side only. The narrow bridge
  # removes the arrow-marked corner; the warm model begins directly at splice_x.
  y_final <- y_reference
  y_final[bridge_idx] <- bridge_y
  
  warm_idx <- x_seq > splice_x
  y_final[warm_idx] <- warm_fit$y[warm_idx]
  
  # Enforce monotonicity only from the bridge start onward. Values colder than
  # bridge_start_x remain bit-for-bit identical to the original reference curve.
  y_final[bridge_start_index:length(y_final)] <- cummax(
    pmax(
      y_final[bridge_start_index:length(y_final)],
      y_final[bridge_start_index]
    )
  )
  
  y_final <- pmin(pmax(y_final, 0), 1)
  
  exact_left_index <- seq_len(max(1L, bridge_start_index - 1L))
  max_left_difference <- max(
    abs(y_final[exact_left_index] - y_reference[exact_left_index]),
    na.rm = TRUE
  )
  
  if (!is.finite(max_left_difference) || max_left_difference > 1e-12) {
    stop(
      "Far cold-side preservation check failed; maximum difference = ",
      format(max_left_difference, scientific = TRUE),
      call. = FALSE
    )
  }
  
  curve_df <- tibble(
    bio1 = x_seq,
    suit01 = y_final,
    suit01_reference = y_reference,
    preserved_left = seq_along(x_seq) < bridge_start_index,
    threshold_bridge = seq_along(x_seq) >= bridge_start_index &
      seq_along(x_seq) <= splice_index
  )
  
  # Recalculate the reported threshold from the FINAL plotted curve. This
  # guarantees that the red line is the strict maximum-positive-curvature point
  # of the curve that is actually shown, rather than a retained coordinate from
  # the earlier isotonic reference curve.
  final_derivatives <- find_coldside_curvature_threshold(curve_df)
  threshold_x <- final_derivatives$threshold
  threshold_index <- final_derivatives$threshold_index
  
  curve_df <- curve_df %>%
    mutate(
      slope = final_derivatives$slope,
      second_derivative = final_derivatives$second_derivative,
      normalized_slope = final_derivatives$normalized_slope,
      normalized_second_derivative =
        final_derivatives$normalized_second_derivative,
      signed_curvature = final_derivatives$signed_curvature,
      threshold_candidate = final_derivatives$candidate_mask,
      is_threshold = ifelse(
        is.na(threshold_index),
        FALSE,
        row_number() == threshold_index
      )
    )
  
  envelope_points$envelope_reference <- as.numeric(
    approx(
      x = x_seq,
      y = y_reference,
      xout = envelope_points$bio1,
      rule = 2
    )$y
  )
  
  envelope_points$envelope_final <- as.numeric(
    approx(
      x = x_seq,
      y = y_final,
      xout = envelope_points$bio1,
      rule = 2
    )$y
  )
  
  envelope_points$preserved_left <-
    envelope_points$bio1 < bridge_start_x
  
  envelope_points$threshold_bridge <-
    envelope_points$bio1 >= bridge_start_x &
    envelope_points$bio1 <= splice_x
  
  list(
    method = paste0(
      "cold density-adaptive envelope; fixed 0.1 °C all-point warm q99; ",
      "exact far-cold baseline; ", bridge_method, "; ",
      warm_fit$method
    ),
    threshold_method = paste0(
      final_derivatives$threshold_method,
      "; computed on final plotted curve"
    ),
    threshold = threshold_x,
    threshold_index = threshold_index,
    threshold_reference = threshold_seed_x,
    threshold_search_stage = final_derivatives$search_stage,
    threshold_curvature = final_derivatives$threshold_curvature,
    maximum_candidate_curvature =
      final_derivatives$maximum_candidate_curvature,
    threshold_is_strict_curvature_max =
      final_derivatives$strict_curvature_max,
    bridge_start_temperature = bridge_start_x,
    splice_temperature = splice_x,
    blend_end_temperature = splice_x,
    bridge_baseline_fraction_used = bridge_baseline_fraction_use,
    bridge_min_left_width_used = bridge_min_left_width_use,
    bridge_max_left_width_used = bridge_max_left_width_use,
    bridge_left_slope_estimated = bridge_left_slope_estimated,
    bridge_left_slope_used = bridge_left_slope,
    bridge_method = bridge_method,
    upper_anchor = upper_anchor,
    max_left_difference = max_left_difference,
    warm_fit_method = warm_fit$method,
    warm_fit_parameters = warm_fit$parameters,
    curve_df = curve_df,
    envelope_points = envelope_points,
    warm_q99_points = warm_q99_points,
    smooth_fit = smooth_fit
  )
}

# ------------------------------------------------------------
# 11. Derivatives and normalized cold-side curvature threshold
# Same logic as the revised Antarctic method
# ------------------------------------------------------------

central_diff <- function(x, y) {
  n <- length(x)
  out <- rep(NA_real_, n)
  
  if (n < 3) return(out)
  
  out[1] <- (y[2] - y[1]) / (x[2] - x[1])
  out[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
  out[2:(n - 1)] <-
    (y[3:n] - y[1:(n - 2)]) /
    (x[3:n] - x[1:(n - 2)])
  
  out
}

find_coldside_curvature_threshold <- function(
    curve_df,
    edge_trim = threshold_edge_trim,
    rise_trim = threshold_rise_trim,
    min_slope_fraction = threshold_min_slope_fraction
) {
  x <- as.numeric(curve_df$bio1)
  y <- as.numeric(curve_df$suit01)
  
  d1_raw <- central_diff(x, y)
  d2_raw <- central_diff(x, d1_raw)
  
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  y_min <- min(y, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  
  x_rng <- x_max - x_min
  y_rng <- y_max - y_min
  
  empty_result <- function(method_text) {
    list(
      threshold = NA_real_,
      threshold_index = NA_integer_,
      slope = d1_raw,
      second_derivative = d2_raw,
      normalized_slope = rep(NA_real_, length(x)),
      normalized_second_derivative = rep(NA_real_, length(x)),
      signed_curvature = rep(NA_real_, length(x)),
      candidate_mask = rep(FALSE, length(x)),
      search_stage = "none",
      threshold_curvature = NA_real_,
      maximum_candidate_curvature = NA_real_,
      strict_curvature_max = FALSE,
      threshold_method = method_text
    )
  }
  
  if (!is.finite(x_rng) || x_rng <= 0 ||
      !is.finite(y_rng) || y_rng <= 0) {
    return(empty_result(
      "strict maximum positive curvature on normalized cold-side envelope"
    ))
  }
  
  x01 <- (x - x_min) / x_rng
  y01 <- (y - y_min) / y_rng
  
  d1_norm <- central_diff(x01, y01)
  d2_norm <- central_diff(x01, d1_norm)
  curvature_norm <- d2_norm / ((1 + d1_norm^2)^(3 / 2))
  
  x_low  <- edge_trim
  x_high <- 1 - edge_trim
  y_low  <- rise_trim[1]
  y_high <- rise_trim[2]
  
  slope_core <- which(
    is.finite(d1_norm) &
      d1_norm > 0 &
      x01 >= x_low &
      x01 <= x_high
  )
  
  if (length(slope_core) == 0) {
    return(empty_result(
      "strict maximum positive curvature on normalized cold-side envelope"
    ))
  }
  
  peak_slope_index <- slope_core[which.max(d1_norm[slope_core])]
  x_peak_slope <- x01[peak_slope_index]
  max_positive_slope <- d1_norm[peak_slope_index]
  slope_floor <- min_slope_fraction * max_positive_slope
  
  # Stage 1: the intended method. Search the accelerating cold-side branch,
  # before the maximum slope, within the lower 0.5%-30% of total response rise.
  cand <- which(
    is.finite(d1_norm) &
      is.finite(d2_norm) &
      is.finite(curvature_norm) &
      d1_norm >= slope_floor &
      d2_norm > 0 &
      curvature_norm > 0 &
      y01 >= y_low &
      y01 <= y_high &
      x01 >= x_low &
      x01 <= min(x_high, x_peak_slope)
  )
  search_stage <- "primary_0.5_to_30_percent_rise"
  
  # Stage 2 keeps the method curvature-based: only relax the response-height
  # restriction; positive curvature and pre-maximum-slope restrictions remain.
  if (length(cand) == 0) {
    cand <- which(
      is.finite(d1_norm) &
        is.finite(d2_norm) &
        is.finite(curvature_norm) &
        d1_norm >= slope_floor &
        d2_norm > 0 &
        curvature_norm > 0 &
        x01 >= x_low &
        x01 <= min(x_high, x_peak_slope)
    )
    search_stage <- "relaxed_height_strict_curvature"
  }
  
  # Stage 3 still selects maximum positive curvature. It relaxes only the slope
  # floor so a very gradual cold-side departure does not return T = NA.
  if (length(cand) == 0) {
    cand <- which(
      is.finite(d1_norm) &
        is.finite(d2_norm) &
        is.finite(curvature_norm) &
        d1_norm > 0 &
        d2_norm > 0 &
        curvature_norm > 0 &
        x01 >= x_low &
        x01 <= min(x_high, x_peak_slope)
    )
    search_stage <- "relaxed_slope_strict_curvature"
  }
  
  candidate_mask <- rep(FALSE, length(x))
  if (length(cand) > 0) candidate_mask[cand] <- TRUE
  
  threshold_index <- if (length(cand) > 0) {
    cand[which.max(curvature_norm[cand])]
  } else {
    NA_integer_
  }
  
  threshold <- if (is.na(threshold_index)) NA_real_ else x[threshold_index]
  threshold_curvature <- if (is.na(threshold_index)) {
    NA_real_
  } else {
    curvature_norm[threshold_index]
  }
  
  maximum_candidate_curvature <- if (length(cand) > 0) {
    max(curvature_norm[cand], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  strict_curvature_max <- is.finite(threshold_curvature) &&
    is.finite(maximum_candidate_curvature) &&
    isTRUE(all.equal(
      threshold_curvature,
      maximum_candidate_curvature,
      tolerance = 1e-12
    ))
  
  list(
    threshold = threshold,
    threshold_index = threshold_index,
    slope = d1_raw,
    second_derivative = d2_raw,
    normalized_slope = d1_norm,
    normalized_second_derivative = d2_norm,
    signed_curvature = curvature_norm,
    candidate_mask = candidate_mask,
    search_stage = search_stage,
    threshold_curvature = threshold_curvature,
    maximum_candidate_curvature = maximum_candidate_curvature,
    strict_curvature_max = strict_curvature_max,
    threshold_method = paste0(
      "strict maximum positive curvature on normalized cold-side ",
      "density-adaptive envelope; before maximum slope; search = ",
      search_stage
    )
  )
}

# ------------------------------------------------------------
# 12. Plot helpers
# Same plot parameters as the Antarctic panels
# ------------------------------------------------------------

theme_response <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
      plot.subtitle = element_text(
        size = 8.5,
        hjust = 0.5,
        colour = "black"
      ),
      axis.title = element_text(size = 11, colour = "black"),
      axis.text = element_text(size = 9.5, colour = "black"),
      axis.line = element_line(linewidth = 0.6, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

make_plot_sample <- function(dat, n = max_scatter_n) {
  if (nrow(dat) <= n) return(dat)
  dat %>% slice_sample(n = n)
}

make_one_panel <- function(
    sample_out,
    fit_obj,
    group_name,
    xlim_use
) {
  scatter_dat <- make_plot_sample(sample_out)
  
  thr <- fit_obj$threshold
  x_range <- diff(xlim_use)
  
  label_x <- if (is.finite(thr)) {
    min(
      thr + 0.035 * x_range,
      xlim_use[2] - 0.12 * x_range
    )
  } else {
    xlim_use[1] + 0.05 * x_range
  }
  
  p <- ggplot() +
    geom_point(
      data = scatter_dat,
      aes(x = bio1, y = suit01_clip),
      colour = point_col,
      size = point_size,
      alpha = point_alpha
    ) +
    geom_line(
      data = fit_obj$curve_df,
      aes(x = bio1, y = suit01),
      colour = curve_col,
      linewidth = curve_linewidth
    )
  
  # Add a threshold line only when a finite threshold is available, avoiding
  # geom_vline() warnings caused by xintercept = NA.
  if (is.finite(thr)) {
    p <- p +
      geom_vline(
        xintercept = thr,
        colour = threshold_col,
        linewidth = threshold_width,
        linetype = threshold_type,
        alpha = threshold_alpha
      )
  }
  
  p +
    annotate(
      "text",
      x = label_x,
      y = 0.08,
      label = if (is.finite(thr)) {
        paste0("T = ", sprintf("%.2f", thr))
      } else {
        "T = NA"
      },
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      colour = "black"
    ) +
    coord_cartesian(
      xlim = xlim_use,
      ylim = c(0, 1.10),
      clip = "off"
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1.0),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = group_name,
      subtitle = "cold density threshold; warm 0.1°C all-point q99 edge",
      x = "Annual mean temperature (°C)",
      y = "Suitability (0–1)"
    ) +
    theme_response()
}

# ------------------------------------------------------------
# 13. Build one Arctic figure: Growth form or Color type
# ------------------------------------------------------------

make_arctic_response_figure <- function(
    group_type_name,
    group_names_use
) {
  info_tbl <- tibble(
    pole = "Arctic",
    group_type = group_type_name,
    group_name = group_names_use,
    suitability_file = vapply(
      group_names_use,
      function(g) make_suitability_path(group_type_name, g),
      character(1)
    ),
    bio1_file = make_bio1_path()
  )
  
  print(info_tbl)
  walk(info_tbl$suitability_file, check_file_exists)
  check_file_exists(info_tbl$bio1_file[1])
  
  res_list <- map(
    info_tbl$group_name,
    function(group_name) {
      message(
        "\nProcessing: Arctic | ",
        group_type_name, " | ", group_name
      )
      
      sample_out <- build_arctic_sample(
        group_type_name = group_type_name,
        group_name = group_name
      )
      
      fit_obj <- fit_density_adaptive_envelope(
        sample_out = sample_out,
        group_type_name = group_type_name,
        group_name = group_name
      )
      
      # Save diagnostic local-envelope points for checking how each local
      # temperature window was classified.
      write.csv(
        fit_obj$envelope_points,
        file.path(
          out_dir,
          paste0(
            "Cold_local_envelope_points_Arctic_",
            gsub(" ", "_", group_type_name), "_",
            group_name,
            ".csv"
          )
        ),
        row.names = FALSE
      )
      
      write.csv(
        fit_obj$warm_q99_points,
        file.path(
          out_dir,
          paste0(
            "Warm_fixed_0.1C_q99_points_Arctic_",
            gsub(" ", "_", group_type_name), "_",
            group_name,
            ".csv"
          )
        ),
        row.names = FALSE
      )
      
      list(
        group_name = group_name,
        sample_out = sample_out,
        fit_obj = fit_obj,
        threshold = fit_obj$threshold
      )
    }
  )
  
  all_x <- unlist(map(res_list, ~ .x$sample_out$bio1))
  xlim_use <- range(all_x, na.rm = TRUE)
  x_pad <- diff(xlim_use) * 0.05
  xlim_use <- c(xlim_use[1] - x_pad, xlim_use[2] + x_pad)
  
  plot_list <- map(
    res_list,
    ~ make_one_panel(
      sample_out = .x$sample_out,
      fit_obj = .x$fit_obj,
      group_name = .x$group_name,
      xlim_use = xlim_use
    )
  )
  
  threshold_tbl <- map_dfr(
    res_list,
    function(z) {
      env <- z$fit_obj$envelope_points
      
      tibble(
        pole = "Arctic",
        group_type = group_type_name,
        group_name = z$group_name,
        threshold_bio1 = z$threshold,
        reference_threshold_bio1 = z$fit_obj$threshold_reference,
        threshold_search_stage = z$fit_obj$threshold_search_stage,
        threshold_curvature = z$fit_obj$threshold_curvature,
        maximum_candidate_curvature =
          z$fit_obj$maximum_candidate_curvature,
        threshold_is_strict_curvature_max =
          z$fit_obj$threshold_is_strict_curvature_max,
        method = z$fit_obj$method,
        threshold_method = z$fit_obj$threshold_method,
        sampling_mode = arctic_sampling_mode,
        sample_n = nrow(z$sample_out),
        cold_local_envelope_point_n = nrow(env),
        warm_fixed_q99_point_n = nrow(z$fit_obj$warm_q99_points),
        clustered_window_fraction = mean(
          env$clustered,
          na.rm = TRUE
        ),
        median_cluster_weight = median(
          env$cluster_weight,
          na.rm = TRUE
        ),
        median_window_n = median(
          env$n_window,
          na.rm = TRUE
        ),
        temp_step = temp_step,
        initial_window_width = 2 * window_half_width_initial,
        maximum_window_width = 2 * window_half_width_max,
        local_quantile = local_quantile,
        cluster_mass = cluster_mass,
        cluster_width_full = cluster_width_full,
        cluster_width_none = cluster_width_none,
        envelope_spar = envelope_spar,
        bridge_start_temperature = z$fit_obj$bridge_start_temperature,
        splice_temperature = z$fit_obj$splice_temperature,
        bridge_baseline_fraction_used =
          z$fit_obj$bridge_baseline_fraction_used,
        bridge_min_left_width_used =
          z$fit_obj$bridge_min_left_width_used,
        bridge_max_left_width_used =
          z$fit_obj$bridge_max_left_width_used,
        bridge_left_slope_estimated =
          z$fit_obj$bridge_left_slope_estimated,
        bridge_left_slope_used =
          z$fit_obj$bridge_left_slope_used,
        bridge_method = z$fit_obj$bridge_method,
        blend_end_temperature = z$fit_obj$blend_end_temperature,
        upper_anchor = z$fit_obj$upper_anchor,
        warm_fit_method = z$fit_obj$warm_fit_method,
        max_left_difference = z$fit_obj$max_left_difference,
        left_preserve_after_threshold = left_preserve_after_threshold,
        left_preserve_rise_fraction = left_preserve_rise_fraction,
        left_preserve_max_after_threshold =
          left_preserve_max_after_threshold,
        threshold_bridge_baseline_fraction =
          threshold_bridge_baseline_fraction,
        threshold_bridge_min_left_width =
          threshold_bridge_min_left_width,
        threshold_bridge_max_left_width =
          threshold_bridge_max_left_width,
        warm_bin_width = warm_bin_width,
        warm_min_bin_n = warm_min_bin_n,
        warm_quantile = warm_quantile,
        warm_start_after_threshold = warm_start_after_threshold,
        warm_underfit_penalty = warm_underfit_penalty,
        warm_spline_spar = warm_spline_spar,
        warm_derivative_spar = warm_derivative_spar,
        warm_fit_upper = unname(z$fit_obj$warm_fit_parameters["upper"]),
        warm_fit_rate = unname(z$fit_obj$warm_fit_parameters["rate"]),
        warm_fit_shape = unname(z$fit_obj$warm_fit_parameters["shape"])
      )
    }
  )
  
  write.csv(
    threshold_tbl,
    file.path(
      out_dir,
      paste0(
        "Thermal_thresholds_Arctic_",
        gsub(" ", "_", group_type_name),
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  p_final <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = paste0(
        "Arctic | ",
        group_type_name,
        " thermal response curves"
      ),
      subtitle = ""
    ) &
    theme(
      plot.title = element_text(
        size = 14,
        hjust = 0.5,
        face = "plain"
      ),
      plot.subtitle = element_text(
        size = 10,
        hjust = 0.5
      )
    )
  
  list(
    plot = p_final,
    thresholds = threshold_tbl,
    results = res_list
  )
}

# ------------------------------------------------------------
# 14. Calculate and plot Arctic groups only
# ------------------------------------------------------------

arc_growth <- make_arctic_response_figure(
  group_type_name = "Growth form",
  group_names_use = growth_groups
)

arc_colour <- make_arctic_response_figure(
  group_type_name = "Color type",
  group_names_use = colour_groups
)

p_arc_growth <- arc_growth$plot
p_arc_colour <- arc_colour$plot

p_arc_growth
p_arc_colour

# ------------------------------------------------------------
# 15. Save figures
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    out_dir,
    "Fig_response_Arctic_Growth_form_cold_density_warm_q99.png"
  ),
  plot = p_arc_growth,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = figure_dpi,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Fig_response_Arctic_Color_type_cold_density_warm_q99.png"
  ),
  plot = p_arc_colour,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = figure_dpi,
  bg = "white"
)

# ------------------------------------------------------------
# 16. Combined Arctic threshold table
# ------------------------------------------------------------

threshold_all_arctic <- bind_rows(
  arc_growth$thresholds,
  arc_colour$thresholds
)

write.csv(
  threshold_all_arctic,
  file.path(
    out_dir,
    "Thermal_thresholds_Arctic_all_groups_cold_density_warm_q99.csv"
  ),
  row.names = FALSE
)

threshold_all_arctic