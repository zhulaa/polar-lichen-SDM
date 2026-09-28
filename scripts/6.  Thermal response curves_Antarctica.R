############################################################
# Antarctic thermal response curves from original current suitability rasters
#
# Independent Antarctic-only script extracted from:
# “响应曲线包络线（南极表现较好）(1).R”
#
# Data source:
# RESULT/Antarctic/<Group type>/<Group>/current.img
# EnvironmentData/Antarctic/current/bio1.asc
#
# Suitability raster range: 0–1000
# suit01 = suitability_raw / 1000
#
# Sampling:
# spatial block + five suitability strata
#
# Fitting:
# direct thermal response fitted to all sampled suitability values;
# the final response is lightly smoothed and constrained to be non-decreasing.
#
# Threshold:
# maximum positive curvature on the normalized cold-side response,
# searched before the maximum-slope point.
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

mask_path <- "E:/1 MY PROJECTS/Climate Change/antarctic_non_permanent_glacier_mask/非永久冰川区.tif"

out_dir <- file.path(
  base_work_dir,
  "_response_curve_figures_Antarctic_direct_coldside_curvature"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Groups and settings
# ------------------------------------------------------------

growth_groups <- c("Crustose", "Foliose", "Fruticose")
colour_groups <- c("Pale", "Bright", "Dark")

suit_breaks <- c(0, 200, 400, 600, 800, 1000)
suit_labels <- c("C1", "C2", "C3", "C4", "C5")

block_size_ant <- 3

# Antarctic direct-response method
direct_df <- 3
direct_bin_width <- 0.10
direct_min_bin_n <- 3
spline_spar_direct <- 0.78

# Curve grid and threshold search
grid_n <- 800
ant_edge_trim <- 0.02
ant_rise_trim <- c(0.005, 0.30)
ant_min_slope_fraction <- 0.03

# Plotting
max_scatter_n <- 5000
curve_col <- "black"
threshold_col <- "#E53935"
point_col <- "grey55"

set.seed(123)

# ------------------------------------------------------------
# 3. Path helpers
# ------------------------------------------------------------

make_suitability_path <- function(group_type_name, group_name) {
  file.path(
    base_result_dir,
    "Antarctic",
    group_type_name,
    group_name,
    "current.img"
  )
}

make_bio1_path <- function() {
  file.path(
    base_env_dir,
    "Antarctic",
    "current",
    "bio1.asc"
  )
}

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File not found:\n", path, call. = FALSE)
  }
}

# ------------------------------------------------------------
# 4. Antarctic raster alignment
# ------------------------------------------------------------

read_project_mask_ant <- function(r_path, mask_path, method = "bilinear") {
  check_file_exists(r_path)
  check_file_exists(mask_path)
  
  r <- terra::rast(r_path)
  m <- terra::rast(mask_path)
  
  r_p <- terra::project(r, m, method = method)
  terra::mask(r_p, m)
}

# ------------------------------------------------------------
# 5. Build original block + suitability-stratified Antarctic sample
# ------------------------------------------------------------

build_current_sample_ant <- function(group_type_name, group_name) {
  su_path <- make_suitability_path(group_type_name, group_name)
  bio1_path <- make_bio1_path()
  
  check_file_exists(su_path)
  check_file_exists(bio1_path)
  
  su_current <- read_project_mask_ant(
    r_path = su_path,
    mask_path = mask_path,
    method = "bilinear"
  )
  
  bio1_r <- read_project_mask_ant(
    r_path = bio1_path,
    mask_path = mask_path,
    method = "bilinear"
  )
  
  names(su_current) <- "suitability_raw"
  names(bio1_r) <- "bio1"
  
  cell_id <- which(!is.na(terra::values(su_current)))
  xy <- terra::xyFromCell(su_current, cell_id)
  val <- terra::values(su_current, mat = FALSE)[cell_id]
  
  su_df <- data.frame(
    cell = cell_id,
    x = xy[, 1],
    y = xy[, 2],
    suitability_raw = as.numeric(val)
  ) %>%
    filter(
      !is.na(suitability_raw),
      is.finite(suitability_raw),
      suitability_raw >= 0,
      suitability_raw <= 1000
    ) %>%
    mutate(
      suit01 = suitability_raw / 1000,
      suit01_clip = pmin(pmax(suit01, 0.001), 0.999),
      su_class5 = cut(
        suitability_raw,
        breaks = suit_breaks,
        include.lowest = TRUE,
        labels = suit_labels
      )
    ) %>%
    filter(!is.na(su_class5))
  
  rc <- terra::rowColFromCell(su_current, su_df$cell)
  su_df$row <- rc[, 1]
  su_df$col <- rc[, 2]
  
  su_df <- su_df %>%
    mutate(
      block_row = floor((row - 1) / block_size_ant),
      block_col = floor((col - 1) / block_size_ant),
      block_id = paste(block_row, block_col, sep = "_")
    )
  
  sample_df <- su_df %>%
    group_by(block_id, su_class5) %>%
    mutate(rand = runif(n())) %>%
    arrange(rand, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup() %>%
    select(-rand)
  
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
  
  # WorldClim bio1 may be stored as °C × 10
  if (median(abs(sample_out$bio1), na.rm = TRUE) > 80) {
    message(
      "Converting bio1 from °C × 10 to °C for Antarctic | ",
      group_type_name, " | ", group_name
    )
    sample_out <- sample_out %>% mutate(bio1 = bio1 / 10)
  }
  
  message(
    "\nPrepared sample: Antarctic | ", group_type_name, " | ", group_name,
    "\n  n = ", nrow(sample_out),
    "\n  suitability_raw range = ",
    paste(round(range(sample_out$suitability_raw, na.rm = TRUE), 2), collapse = " to "),
    "\n  suit01 range = ",
    paste(round(range(sample_out$suit01, na.rm = TRUE), 3), collapse = " to "),
    "\n  bio1 range = ",
    paste(round(range(sample_out$bio1, na.rm = TRUE), 3), collapse = " to ")
  )
  
  sample_out
}

# ------------------------------------------------------------
# 6. Antarctic direct-response fitting and threshold helpers
# ------------------------------------------------------------

safe_direct_fit <- function(formula_obj, dat, weight_col = NULL) {
  
  dat2 <- dat
  
  if (!is.null(weight_col) && weight_col %in% names(dat2)) {
    dat2$.__w__ <- as.numeric(dat2[[weight_col]])
  } else {
    dat2$.__w__ <- 1
  }
  
  dat2 <- dat2 %>%
    dplyr::filter(
      !is.na(.__w__),
      is.finite(.__w__),
      .__w__ > 0
    )
  
  fit_glm <- tryCatch(
    glm(
      formula_obj,
      data = dat2,
      weights = .__w__,
      family = quasibinomial(link = "logit")
    ),
    error = function(e) {
      message("GLM failed: ", conditionMessage(e))
      NULL
    }
  )
  
  if (!is.null(fit_glm)) return(fit_glm)
  
  fit_lm <- tryCatch(
    lm(
      formula_obj,
      data = dat2,
      weights = .__w__
    ),
    error = function(e) {
      message("LM failed: ", conditionMessage(e))
      NULL
    }
  )
  
  fit_lm
}

make_ns_term_local <- function(var, dat, df_target = 3) {
  
  x <- dat[[var]]
  x <- x[!is.na(x) & is.finite(x)]
  nu <- length(unique(round(x, 6)))
  
  if (nu <= 4) {
    return(var)
  }
  
  use_df <- min(df_target, nu - 1)
  paste0("splines::ns(", var, ", df = ", use_df, ")")
}

safe_predict_curve <- function(fit, newdata) {
  
  pred <- tryCatch(
    {
      if (inherits(fit, "glm")) {
        as.numeric(predict(fit, newdata = newdata, type = "response"))
      } else {
        as.numeric(predict(fit, newdata = newdata))
      }
    },
    error = function(e) rep(NA_real_, nrow(newdata))
  )
  
  pmin(pmax(pred, 0), 1)
}

central_diff <- function(x, y) {
  
  n <- length(x)
  out <- rep(NA_real_, n)
  
  if (n < 3) return(out)
  
  out[1] <- (y[2] - y[1]) / (x[2] - x[1])
  out[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
  out[2:(n - 1)] <- (y[3:n] - y[1:(n - 2)]) /
    (x[3:n] - x[1:(n - 2)])
  
  out
}

find_threshold_ant_curvature <- function(
    curve_df,
    edge_trim = ant_edge_trim,
    rise_trim = ant_rise_trim,
    min_slope_fraction = ant_min_slope_fraction
) {
  
  x <- as.numeric(curve_df$bio1)
  y <- as.numeric(curve_df$suit01)
  
  # Raw derivatives are retained for diagnostics.
  d1_raw <- central_diff(x, y)
  d2_raw <- central_diff(x, d1_raw)
  
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  x_rng <- x_max - x_min
  
  y_min <- min(y, na.rm = TRUE)
  y_max <- max(y, na.rm = TRUE)
  y_rng <- y_max - y_min
  
  if (!is.finite(x_rng) || x_rng <= 0 ||
      !is.finite(y_rng) || y_rng <= 0) {
    return(list(
      threshold = NA_real_,
      threshold_index = NA_integer_,
      slope = d1_raw,
      second_derivative = d2_raw,
      normalized_slope = rep(NA_real_, length(x)),
      normalized_second_derivative = rep(NA_real_, length(x)),
      signed_curvature = rep(NA_real_, length(x)),
      threshold_method = paste0(
        "maximum positive curvature on normalized cold-side direct response; ",
        "before maximum slope"
      )
    ))
  }
  
  # Normalize both axes before calculating curvature. Without this step,
  # temperature spans tens of degrees while suitability spans only 0–1, so the
  # raw-coordinate curvature tends to select a later point near the middle of
  # the sigmoid. The normalized curvature instead captures the early cold-side
  # bend where the response visibly departs from the x-axis.
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
  
  if (length(slope_core) > 0) {
    peak_slope_index <- slope_core[which.max(d1_norm[slope_core])]
    x_peak_slope <- x01[peak_slope_index]
    max_positive_slope <- max(d1_norm[slope_core], na.rm = TRUE)
  } else {
    peak_slope_index <- NA_integer_
    x_peak_slope <- x_high
    max_positive_slope <- NA_real_
  }
  
  slope_floor <- if (is.finite(max_positive_slope)) {
    min_slope_fraction * max_positive_slope
  } else {
    0
  }
  
  # Restrict the search to the lower, accelerating branch of the response and
  # to temperatures colder than the maximum-slope point. This prevents the
  # later central part of the sigmoid from being selected.
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
  
  # First fallback: retain the cold-side and pre-maximum-slope restrictions,
  # but remove the suitability-window restriction.
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
  }
  
  # Final fallback for unusually flat curves.
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
  }
  
  threshold_index <- if (length(cand) > 0) {
    cand[which.max(curvature_norm[cand])]
  } else {
    NA_integer_
  }
  
  threshold <- if (is.na(threshold_index)) NA_real_ else x[threshold_index]
  
  list(
    threshold = threshold,
    threshold_index = threshold_index,
    slope = d1_raw,
    second_derivative = d2_raw,
    normalized_slope = d1_norm,
    normalized_second_derivative = d2_norm,
    signed_curvature = curvature_norm,
    threshold_method = paste0(
      "maximum positive curvature on normalized cold-side direct response; ",
      "before maximum slope"
    )
  )
}

fit_direct_response_curve <- function(
    sample_out,
    df_target = direct_df,
    bin_width = direct_bin_width,
    min_bin_n = direct_min_bin_n
) {
  
  dat <- sample_out %>%
    transmute(
      bio1 = as.numeric(bio1),
      suit01_clip = as.numeric(suit01_clip)
    ) %>%
    filter(
      !is.na(bio1),
      !is.na(suit01_clip),
      is.finite(bio1),
      is.finite(suit01_clip)
    ) %>%
    mutate(
      suit01_clip = pmin(pmax(suit01_clip, 0.001), 0.999)
    )
  
  if (nrow(dat) < 10) {
    stop("Too few data points for Antarctic direct response fitting.")
  }
  
  thermal_term <- make_ns_term_local("bio1", dat, df_target = df_target)
  fml <- as.formula(paste("suit01_clip ~", thermal_term))
  fit <- safe_direct_fit(fml, dat)
  
  x_seq <- seq(
    min(dat$bio1, na.rm = TRUE),
    max(dat$bio1, na.rm = TRUE),
    length.out = grid_n
  )
  
  if (!is.null(fit)) {
    curve_df <- tibble(
      bio1 = x_seq,
      suit01 = safe_predict_curve(
        fit,
        newdata = data.frame(bio1 = x_seq)
      )
    )
    method_use <- "direct thermal response fit"
  } else {
    message("Direct model failed; falling back to binned direct smoothing.")
    
    fit_pts <- dat %>%
      mutate(temp_bin = round(bio1 / bin_width) * bin_width) %>%
      group_by(temp_bin) %>%
      filter(n() >= min_bin_n) %>%
      summarise(
        bio1 = mean(bio1, na.rm = TRUE),
        suit01_clip = mean(suit01_clip, na.rm = TRUE),
        n_all = n(),
        .groups = "drop"
      ) %>%
      arrange(bio1)
    
    if (nrow(fit_pts) < 6) {
      stop("Both direct model and binned smoothing failed: too few usable temperature bins.")
    }
    
    sp <- smooth.spline(
      x = fit_pts$bio1,
      y = fit_pts$suit01_clip,
      w = sqrt(fit_pts$n_all),
      spar = spline_spar_direct
    )
    
    curve_df <- tibble(
      bio1 = x_seq,
      suit01 = predict(sp, x = x_seq)$y
    )
    method_use <- "binned direct thermal response fit"
  }
  
  # Retain the original final light smoothing and cold-side monotonic rule.
  sp2 <- tryCatch(
    smooth.spline(
      x = curve_df$bio1,
      y = curve_df$suit01,
      spar = 0.45
    ),
    error = function(e) NULL
  )
  
  if (!is.null(sp2)) {
    curve_df$suit01 <- predict(sp2, x = curve_df$bio1)$y
  }
  
  curve_df <- curve_df %>%
    mutate(suit01 = pmin(pmax(suit01, 0), 1))
  
  curve_df$suit01 <- cummax(curve_df$suit01)
  curve_df$suit01 <- pmin(pmax(curve_df$suit01, 0), 1)
  
  thr_obj <- find_threshold_ant_curvature(
    curve_df = curve_df,
    edge_trim = ant_edge_trim,
    rise_trim = ant_rise_trim,
    min_slope_fraction = ant_min_slope_fraction
  )
  
  curve_df <- curve_df %>%
    mutate(
      slope = thr_obj$slope,
      second_derivative = thr_obj$second_derivative,
      normalized_slope = thr_obj$normalized_slope,
      normalized_second_derivative = thr_obj$normalized_second_derivative,
      signed_curvature = thr_obj$signed_curvature,
      is_threshold = row_number() == thr_obj$threshold_index
    )
  
  list(
    method = method_use,
    fit = fit,
    fit_points = dat,
    curve_df = curve_df,
    threshold = thr_obj$threshold,
    threshold_index = thr_obj$threshold_index,
    threshold_method = thr_obj$threshold_method
  )
}

# ------------------------------------------------------------
# 7. Plot helpers
# ------------------------------------------------------------

theme_response <- function() {
  theme_classic(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(size = 11, face = "plain", hjust = 0.5),
      plot.subtitle = element_text(size = 8.5, hjust = 0.5, colour = "black"),
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

make_one_panel_ant <- function(sample_out, fit_obj, group_name, xlim_use) {
  scatter_dat <- make_plot_sample(sample_out)
  thr <- fit_obj$threshold
  x_range <- diff(xlim_use)
  
  label_x <- ifelse(
    is.na(thr),
    xlim_use[1] + 0.05 * x_range,
    min(thr + 0.035 * x_range, xlim_use[2] - 0.12 * x_range)
  )
  
  p <- ggplot() +
    geom_point(
      data = scatter_dat,
      aes(x = bio1, y = suit01_clip),
      colour = point_col,
      size = 0.85,
      alpha = 0.32
    ) +
    geom_line(
      data = fit_obj$curve_df,
      aes(x = bio1, y = suit01),
      colour = curve_col,
      linewidth = 1.10
    )
  
  if (is.finite(thr)) {
    p <- p +
      geom_vline(
        xintercept = thr,
        colour = threshold_col,
        linewidth = 0.85,
        linetype = "longdash",
        alpha = 0.95
      )
  }
  
  p +
    annotate(
      "text",
      x = label_x,
      y = 0.08,
      label = ifelse(is.na(thr), "T = NA", paste0("T = ", sprintf("%.2f", thr))),
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
      subtitle = paste0(
        "direct smoothed response; threshold = ",
        "cold-side maximum positive curvature"
      ),
      x = "Annual mean temperature (°C)",
      y = "Suitability (0–1)"
    ) +
    theme_response()
}

# ------------------------------------------------------------
# 8. Build one Antarctic figure
# ------------------------------------------------------------

make_response_figure_ant <- function(group_type_name, group_names_use) {
  info_tbl <- tibble(
    pole = "Antarctic",
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
  walk(info_tbl$bio1_file, check_file_exists)
  
  res_list <- map(
    info_tbl$group_name,
    function(group_name) {
      message("\nProcessing: Antarctic | ", group_type_name, " | ", group_name)
      
      sample_out <- build_current_sample_ant(
        group_type_name = group_type_name,
        group_name = group_name
      )
      
      fit_obj <- fit_direct_response_curve(sample_out)
      
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
    ~ make_one_panel_ant(
      sample_out = .x$sample_out,
      fit_obj = .x$fit_obj,
      group_name = .x$group_name,
      xlim_use = xlim_use
    )
  )
  
  threshold_tbl <- map_dfr(
    res_list,
    ~ tibble(
      pole = "Antarctic",
      group_type = group_type_name,
      group_name = .x$group_name,
      threshold_bio1 = .x$threshold,
      method = .x$fit_obj$method,
      threshold_method = .x$fit_obj$threshold_method,
      sample_n = nrow(.x$sample_out),
      fitting_n = nrow(.x$fit_obj$fit_points),
      direct_df = direct_df,
      direct_bin_width = direct_bin_width,
      direct_min_bin_n = direct_min_bin_n,
      spline_spar_direct = spline_spar_direct,
      ant_edge_trim = ant_edge_trim,
      ant_rise_lower = ant_rise_trim[1],
      ant_rise_upper = ant_rise_trim[2],
      ant_min_slope_fraction = ant_min_slope_fraction
    )
  )
  
  write.csv(
    threshold_tbl,
    file.path(
      out_dir,
      paste0(
        "Thermal_thresholds_Antarctic_",
        gsub(" ", "_", group_type_name),
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  p_final <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = paste0(
        "Antarctic | ",
        group_type_name,
        " thermal response curves"
      ),
      subtitle = ""
    ) &
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "plain"),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  
  list(
    plot = p_final,
    thresholds = threshold_tbl,
    results = res_list
  )
}

# ------------------------------------------------------------
# 9. Generate Antarctic figures
# ------------------------------------------------------------

ant_growth <- make_response_figure_ant(
  group_type_name = "Growth form",
  group_names_use = growth_groups
)

ant_colour <- make_response_figure_ant(
  group_type_name = "Color type",
  group_names_use = colour_groups
)

p_ant_growth <- ant_growth$plot
p_ant_colour <- ant_colour$plot

p_ant_growth
p_ant_colour

# ------------------------------------------------------------
# 10. Save Antarctic figures
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    out_dir,
    "Fig_response_Antarctic_Growth_form_direct_coldside_curvature.png"
  ),
  plot = p_ant_growth,
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    out_dir,
    "Fig_response_Antarctic_Color_type_direct_coldside_curvature.png"
  ),
  plot = p_ant_colour,
  width = 180,
  height = 85,
  units = "mm",
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# 11. Combined Antarctic threshold table
# ------------------------------------------------------------

threshold_all_ant <- bind_rows(
  ant_growth$thresholds,
  ant_colour$thresholds
)

write.csv(
  threshold_all_ant,
  file.path(
    out_dir,
    "Thermal_thresholds_Antarctic_all_groups.csv"
  ),
  row.names = FALSE
)

threshold_all_ant