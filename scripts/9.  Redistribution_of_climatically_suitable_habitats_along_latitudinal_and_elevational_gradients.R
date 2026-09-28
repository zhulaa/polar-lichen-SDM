# ============================================================
# Figure 6 | Latitude and elevation redistribution of polar lichens
#
# Purpose
# 1) Calculate suitability-weighted absolute latitude and elevation
#    under current and future climates for all 12 functional groups.
# 2) Calculate 5%, 50%, and 95% suitability-weighted quantiles.
# 3) Calculate mean elevation of high-suitability core habitat (P >= 0.8).
# 4) Summarize latitude/elevation shifts for all 2050s/2090s × 4 SSPs.
# 5) For 2090s SSP5-8.5, create weighted distribution profiles along
#    absolute latitude and elevation.
# 6) For 2090s SSP5-8.5, create exact spatial contribution maps showing
#    which pixels contribute to poleward/equatorward and upslope/downslope
#    shifts in the suitability-weighted mean.
#
# Important
# - Uses glacier-masked suitability rasters.
# - Absolute latitude |lat| is used for BOTH poles:
#     positive shift = poleward
#     negative shift = equatorward
# - Elevation:
#     positive shift = upslope
#     negative shift = downslope
#
# Exact weighted mean:
#   Zbar = sum(P_i * Z_i * A_i) / sum(P_i * A_i)
#
# Exact per-cell contribution to the future-current shift:
#   C_i = Z_i * [ (P_f,i A_i / W_f) - (P_c,i A_i / W_c) ]
# where W = sum(P_i A_i).
# Therefore sum(C_i) = Zbar_future - Zbar_current.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(terra)
library(sf)
library(dplyr)

if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
  stop("Package 'rnaturalearth' is required. Install it with install.packages('rnaturalearth').")
}
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)


# ============================================================
# 1. Paths and global settings
# ============================================================

masked_root <- "F:/1 MY PROJECTS/Climate Change/RESULT_glacier_masked"
work_root   <- "F:/1 MY PROJECTS/Climate Change/分析过程"

# If your DEM is still on E:, only change this line.
dem_path <- "F:/1 MY PROJECTS/Climate Change/wc2.1_10m_elev.tif"

study_area_files <- c(
  "Antarctic" = "F:/1 MY PROJECTS/Climate Change/Study area/antarctica.shp",
  "Arctic"    = "F:/1 MY PROJECTS/Climate Change/Study area/Arctic.shp"
)

fig6_root <- file.path(work_root, "_Fig6_latitude_elevation_redistribution")
table_dir <- file.path(fig6_root, "tables")
profile_dir <- file.path(fig6_root, "weighted_profiles")
contribution_raster_dir <- file.path(fig6_root, "shift_contribution_rasters")
contribution_figure_dir <- file.path(fig6_root, "shift_contribution_maps")
summary_figure_dir <- file.path(fig6_root, "summary_figure")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(profile_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contribution_raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contribution_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_figure_dir, recursive = TRUE, showWarnings = FALSE)

probability_divisor <- 1000
alignment_method <- "bilinear"

# High-suitability core threshold
core_threshold <- 0.8

# Weighted-profile bin widths
latitude_bin_width_deg <- 0.5
elevation_bin_width_m <- 100

# Spatial plotting
map_crs <- c(
  "Antarctic" = "EPSG:3031",
  "Arctic"    = "EPSG:3413"
)

map_max_plot_cells <- 180000
map_colour_quantile <- 0.995

# Figure settings
base_font_pt <- 7
base_family <- "Arial"
map_font_family <- "sans"

summary_width_mm <- 180
summary_height_mm <- 155

profile_width_mm <- 180
profile_height_mm <- 205

contribution_width_mm <- 180
contribution_height_mm <- 95

terraOptions(progress = 1, memfrac = 0.75)


# ============================================================
# 2. Functional groups and future scenarios
# ============================================================

groups <- tribble(
  ~pole,       ~group_type,   ~group_name,
  "Antarctic", "Color type",  "Pale",
  "Antarctic", "Color type",  "Bright",
  "Antarctic", "Color type",  "Dark",
  "Antarctic", "Growth form", "Crustose",
  "Antarctic", "Growth form", "Foliose",
  "Antarctic", "Growth form", "Fruticose",
  "Arctic",    "Color type",  "Pale",
  "Arctic",    "Color type",  "Bright",
  "Arctic",    "Color type",  "Dark",
  "Arctic",    "Growth form", "Crustose",
  "Arctic",    "Growth form", "Foliose",
  "Arctic",    "Growth form", "Fruticose"
) %>%
  mutate(
    group_dir = file.path(masked_root, pole, group_type, group_name),
    output_index = sprintf("%02d", row_number())
  )

future_design <- tribble(
  ~period,  ~scenario, ~scenario_label, ~future_name,
  "2050s",  "SSP126",  "SSP1-2.6",      "2050_SSP126.img",
  "2050s",  "SSP245",  "SSP2-4.5",      "2050_SSP245.img",
  "2050s",  "SSP370",  "SSP3-7.0",      "2050_SSP370.img",
  "2050s",  "SSP585",  "SSP5-8.5",      "2050_SSP585.img",
  "2090s",  "SSP126",  "SSP1-2.6",      "2090_SSP126.img",
  "2090s",  "SSP245",  "SSP2-4.5",      "2090_SSP245.img",
  "2090s",  "SSP370",  "SSP3-7.0",      "2090_SSP370.img",
  "2090s",  "SSP585",  "SSP5-8.5",      "2090_SSP585.img"
)

target_period <- "2090s"
target_scenario <- "SSP585"
target_scenario_label <- "SSP5-8.5"


# ============================================================
# 3. Input checks
# ============================================================

if (!file.exists(dem_path)) {
  stop("DEM not found:\n", dem_path)
}

if (any(!file.exists(study_area_files))) {
  stop(
    "Missing study-area shapefile(s):\n",
    paste(study_area_files[!file.exists(study_area_files)], collapse = "\n")
  )
}

required_files <- c()

for (i in seq_len(nrow(groups))) {
  required_files <- c(
    required_files,
    file.path(groups$group_dir[i], "current.img"),
    file.path(groups$group_dir[i], future_design$future_name)
  )
}

missing_files <- unique(required_files[!file.exists(required_files)])

if (length(missing_files) > 0) {
  stop(
    "Missing suitability raster(s):\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}


# ============================================================
# 4. General helper functions
# ============================================================

raster_range <- function(r) {
  x <- global(r, c("min", "max"), na.rm = TRUE)
  c(
    min = as.numeric(x[1, "min"]),
    max = as.numeric(x[1, "max"])
  )
}


to_prob01 <- function(r, raster_name = "raster") {
  
  rng <- raster_range(r)
  
  if (!all(is.finite(rng))) {
    stop("Raster contains no finite values: ", raster_name)
  }
  
  out <- if (rng["max"] > 1.5) {
    r / probability_divisor
  } else {
    r
  }
  
  rng2 <- raster_range(out)
  
  if (
    rng2["min"] < -1e-6 ||
    rng2["max"] > 1 + 1e-6
  ) {
    stop(
      "Suitability outside 0-1 after rescaling:\n",
      raster_name,
      "\nRange = ",
      rng2["min"],
      " to ",
      rng2["max"]
    )
  }
  
  clamp(out, lower = 0, upper = 1, values = TRUE)
}


align_to_template <- function(
    template_r,
    input_r,
    method = alignment_method
) {
  
  if (isTRUE(compareGeom(template_r, input_r, stopOnError = FALSE))) {
    return(input_r)
  }
  
  if (isTRUE(same.crs(template_r, input_r))) {
    return(
      resample(
        input_r,
        template_r,
        method = method
      )
    )
  }
  
  project(
    input_r,
    template_r,
    method = method,
    mask = FALSE
  )
}


g_sum <- function(r) {
  
  out <- global(
    r,
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  if (
    length(out) == 0 ||
    is.na(out) ||
    !is.finite(out)
  ) {
    return(0)
  }
  
  as.numeric(out)
}


weighted_mean_raster <- function(
    prob_r,
    covar_r,
    area_r
) {
  
  num <- g_sum(
    prob_r * covar_r * area_r
  )
  
  den <- g_sum(
    prob_r * area_r
  )
  
  if (den <= 0) return(NA_real_)
  
  num / den
}


weighted_quantile <- function(
    x,
    w,
    probs = c(0.05, 0.50, 0.95)
) {
  
  ok <- is.finite(x) &
    is.finite(w) &
    w > 0
  
  x <- x[ok]
  w <- w[ok]
  
  if (
    length(x) == 0 ||
    sum(w) <= 0
  ) {
    return(
      setNames(
        rep(NA_real_, length(probs)),
        paste0("q", probs * 100)
      )
    )
  }
  
  ord <- order(x)
  
  x <- x[ord]
  w <- w[ord]
  
  cw <- cumsum(w) / sum(w)
  
  out <- sapply(
    probs,
    function(p) {
      x[which(cw >= p)[1]]
    }
  )
  
  names(out) <- paste0("q", probs * 100)
  
  out
}


weighted_quantile_raster <- function(
    prob_r,
    covar_r,
    area_r,
    probs = c(0.05, 0.50, 0.95)
) {
  
  p <- values(prob_r, mat = FALSE)
  z <- values(covar_r, mat = FALSE)
  a <- values(area_r, mat = FALSE)
  
  w <- p * a
  
  weighted_quantile(
    x = z,
    w = w,
    probs = probs
  )
}


# Robust absolute-latitude raster.
# Works whether the suitability raster is lon/lat or projected.
make_abs_lat_raster <- function(r) {
  
  xy <- crds(
    r,
    df = TRUE,
    na.rm = FALSE
  )
  
  if (is.lonlat(r)) {
    
    abs_lat <- abs(xy[, 2])
    
  } else {
    
    if (is.na(crs(r)) || crs(r) == "") {
      stop("Raster has no CRS; cannot calculate latitude.")
    }
    
    ll <- terra::project(
      as.matrix(xy),
      from = crs(r),
      to = "EPSG:4326"
    )
    
    abs_lat <- abs(ll[, 2])
  }
  
  lat_r <- rast(r)
  values(lat_r) <- abs_lat
  names(lat_r) <- "abs_lat"
  
  lat_r
}


core_mean_elevation <- function(
    prob_r,
    elev_r,
    threshold = core_threshold
) {
  
  core <- ifel(
    prob_r >= threshold,
    1,
    NA
  )
  
  out <- global(
    mask(elev_r, core),
    "mean",
    na.rm = TRUE
  )[1, 1]
  
  if (
    length(out) == 0 ||
    is.na(out) ||
    !is.finite(out)
  ) {
    return(NA_real_)
  }
  
  as.numeric(out)
}


# ============================================================
# 5. Weighted profiles along latitude and elevation
# ============================================================

make_weighted_profile <- function(
    prob_r,
    covar_r,
    area_r,
    bin_width,
    variable_name,
    state_name,
    pole,
    group_type,
    group_name
) {
  
  p <- values(prob_r, mat = FALSE)
  z <- values(covar_r, mat = FALSE)
  a <- values(area_r, mat = FALSE)
  
  w <- p * a
  
  ok <- is.finite(p) &
    is.finite(z) &
    is.finite(a) &
    is.finite(w) &
    w > 0
  
  p <- p[ok]
  z <- z[ok]
  w <- w[ok]
  
  if (length(z) == 0) {
    return(tibble())
  }
  
  bin <- floor(z / bin_width) * bin_width + bin_width / 2
  
  tibble(
    bin = bin,
    weight = w
  ) %>%
    group_by(bin) %>%
    summarise(
      weighted_suitability = sum(weight),
      .groups = "drop"
    ) %>%
    mutate(
      weighted_fraction = weighted_suitability /
        sum(weighted_suitability),
      variable = variable_name,
      state = state_name,
      pole = pole,
      group_type = group_type,
      group_name = group_name
    )
}


# ============================================================
# 6. Exact spatial contribution to weighted-mean shift
# ============================================================

make_shift_contribution_raster <- function(
    current_p,
    future_p,
    covar_r,
    area_r
) {
  
  W_current <- g_sum(
    current_p * area_r
  )
  
  W_future <- g_sum(
    future_p * area_r
  )
  
  if (
    W_current <= 0 ||
    W_future <= 0
  ) {
    stop("Zero WSA denominator while calculating contribution raster.")
  }
  
  contribution <- covar_r * (
    (future_p * area_r / W_future) -
      (current_p * area_r / W_current)
  )
  
  contribution
}


write_float_raster <- function(r, filename) {
  
  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  writeRaster(
    r,
    filename,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      NAflag = -9999,
      gdal = c(
        "COMPRESS=LZW",
        "TILED=YES"
      )
    )
  )
  
  normalizePath(
    filename,
    winslash = "/",
    mustWork = TRUE
  )
}


# ============================================================
# 7. Study-area outlines and Arctic geographic reference
# ============================================================

sf_lines_to_df <- function(sf_obj) {
  
  if (nrow(sf_obj) == 0) {
    return(
      data.frame(
        X = numeric(0),
        Y = numeric(0),
        group_id = character(0)
      )
    )
  }
  
  coords <- as.data.frame(
    st_coordinates(sf_obj)
  )
  
  level_cols <- grep(
    "^L[0-9]+$",
    names(coords),
    value = TRUE
  )
  
  if (length(level_cols) == 0) {
    
    coords$group_id <- "1"
    
  } else {
    
    coords$group_id <- interaction(
      coords[, level_cols, drop = FALSE],
      drop = TRUE,
      lex.order = TRUE
    )
  }
  
  coords
}


study_area_sf <- lapply(
  names(study_area_files),
  function(pole_name) {
    
    x <- st_read(
      study_area_files[[pole_name]],
      quiet = TRUE
    )
    
    x <- st_make_valid(x)
    
    st_transform(
      x,
      map_crs[[pole_name]]
    )
  }
)

names(study_area_sf) <- names(study_area_files)

outline_df_list <- lapply(
  study_area_sf,
  function(x) {
    sf_lines_to_df(
      st_boundary(x)
    )
  }
)


# ------------------------------------------------------------
# Arctic coastlines and international land boundaries
# ------------------------------------------------------------
#
# This follows the Fig. 5 solution:
# - coastline provides complete outlines of North America, Greenland,
#   Europe, Russia and Arctic islands;
# - admin_0_boundary_lines_land provides international land borders.
#
# Natural Earth lines are used directly as sf layers; they are NOT
# converted to polygon boundaries or to a local country polygon file.

cat("
Loading Natural Earth Arctic coastlines and country boundaries...
")

arctic_coastline_sf <- rnaturalearth::ne_download(
  scale = 10,
  type = "coastline",
  category = "physical",
  returnclass = "sf"
)

arctic_admin0_sf <- rnaturalearth::ne_download(
  scale = 10,
  type = "admin_0_boundary_lines_land",
  category = "cultural",
  returnclass = "sf"
)

prepare_arctic_reference <- function(x) {
  
  x <- st_make_valid(x)
  
  # Repair +/-180-degree discontinuities before polar projection.
  if (st_is_longlat(x)) {
    x <- suppressWarnings(
      st_wrap_dateline(
        x,
        options = c(
          "WRAPDATELINE=YES",
          "DATELINEOFFSET=180"
        ),
        quiet = TRUE
      )
    )
  }
  
  x <- st_transform(
    x,
    map_crs[["Arctic"]]
  )
  
  arctic_bbox <- st_bbox(
    study_area_sf[["Arctic"]]
  )
  
  x_pad <- 0.02 * (
    arctic_bbox["xmax"] -
      arctic_bbox["xmin"]
  )
  
  y_pad <- 0.02 * (
    arctic_bbox["ymax"] -
      arctic_bbox["ymin"]
  )
  
  arctic_bbox_expanded <- arctic_bbox
  
  arctic_bbox_expanded["xmin"] <- arctic_bbox["xmin"] - x_pad
  arctic_bbox_expanded["xmax"] <- arctic_bbox["xmax"] + x_pad
  arctic_bbox_expanded["ymin"] <- arctic_bbox["ymin"] - y_pad
  arctic_bbox_expanded["ymax"] <- arctic_bbox["ymax"] + y_pad
  
  x <- suppressWarnings(
    st_crop(
      x,
      arctic_bbox_expanded
    )
  )
  
  x[
    !st_is_empty(x),
  ]
}

arctic_coastline_sf <- prepare_arctic_reference(
  arctic_coastline_sf
)

arctic_admin0_sf <- prepare_arctic_reference(
  arctic_admin0_sf
)

cat(
  "  Coastline features: ",
  nrow(arctic_coastline_sf),
  "
",
sep = ""
)

cat(
  "  International-boundary features: ",
  nrow(arctic_admin0_sf),
  "
",
sep = ""
)


# ============================================================
# 8. Continuous-raster projection for contribution maps
# ============================================================

project_continuous_df_for_plot <- function(
    r,
    pole,
    max_cells = map_max_plot_cells
) {
  
  r2 <- project(
    r,
    map_crs[[pole]],
    method = "bilinear",
    mask = FALSE
  )
  
  if (ncell(r2) > max_cells) {
    
    fact <- ceiling(
      sqrt(
        ncell(r2) / max_cells
      )
    )
    
    r2 <- aggregate(
      r2,
      fact = fact,
      fun = "mean",
      na.rm = TRUE
    )
  }
  
  names(r2) <- "value"
  
  as.data.frame(
    r2,
    xy = TRUE,
    na.rm = TRUE
  )
}


safe_symmetric_limit <- function(
    x,
    probability = map_colour_quantile,
    fallback = 1e-6
) {
  
  x <- abs(
    x[is.finite(x)]
  )
  
  x <- x[x > 0]
  
  if (length(x) == 0) return(fallback)
  
  q <- as.numeric(
    quantile(
      x,
      probs = probability,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  if (!is.finite(q) || q <= 0) fallback else q
}


make_contribution_map <- function(
    raster_file,
    pole,
    outline_df,
    panel_title,
    legend_title
) {
  
  r <- rast(raster_file)
  
  df <- project_continuous_df_for_plot(
    r,
    pole
  )
  
  lim <- safe_symmetric_limit(
    df$value
  )
  
  p <- ggplot(
    df,
    aes(
      x = x,
      y = y,
      fill = value
    )
  ) +
    geom_raster()
  
  # Arctic geographic reference:
  # coastline is lightest; international boundaries slightly darker.
  if (pole == "Arctic") {
    
    p <- p +
      geom_sf(
        data = arctic_coastline_sf,
        inherit.aes = FALSE,
        colour = "grey88",
        linewidth = 0.04
      ) +
      geom_sf(
        data = arctic_admin0_sf,
        inherit.aes = FALSE,
        colour = "grey82",
        linewidth = 0.045
      )
  }
  
  p +
    geom_path(
      data = outline_df,
      aes(
        x = X,
        y = Y,
        group = group_id
      ),
      inherit.aes = FALSE,
      colour = "grey70",
      linewidth = 0.05
    ) +
    coord_sf(
      crs = st_crs(
        map_crs[[pole]]
      ),
      datum = NA,
      expand = FALSE
    ) +
    scale_fill_gradient2(
      low = "#B2182B",
      mid = "white",
      high = "#2166AC",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = squish,
      name = legend_title
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = NULL
    ) +
    theme_void(
      base_size = base_font_pt,
      base_family = map_font_family
    ) +
    theme(
      plot.title = element_text(
        size = base_font_pt,
        face = "bold",
        hjust = 0.5
      ),
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.10
      ),
      legend.position = "bottom",
      legend.title = element_text(
        size = base_font_pt
      ),
      legend.text = element_text(
        size = base_font_pt
      ),
      legend.key.width = unit(
        1.0,
        "cm"
      ),
      plot.margin = margin(
        1,
        1,
        1,
        1
      )
    )
}


# ============================================================
# 9. Main calculation loop
# ============================================================

result_list <- list()
profile_list <- list()
contribution_manifest_list <- list()

result_id <- 1L
profile_id <- 1L
manifest_id <- 1L


for (i in seq_len(nrow(groups))) {
  
  pole_i <- groups$pole[i]
  group_type_i <- groups$group_type[i]
  group_name_i <- groups$group_name[i]
  group_dir_i <- groups$group_dir[i]
  output_index_i <- groups$output_index[i]
  
  cat(
    "\n==================================================\n",
    "Processing: ",
    pole_i,
    " | ",
    group_type_i,
    " | ",
    group_name_i,
    "\n",
    sep = ""
  )
  
  current_file <- file.path(
    group_dir_i,
    "current.img"
  )
  
  current_raw <- rast(
    current_file
  )
  
  current_p_full <- to_prob01(
    current_raw,
    current_file
  )
  
  # DEM aligned to current suitability grid.
  dem_raw <- rast(
    dem_path
  )
  
  dem_current <- align_to_template(
    current_raw,
    dem_raw,
    method = "bilinear"
  )
  
  # Robust absolute latitude for this exact grid.
  abs_lat_current <- make_abs_lat_raster(
    current_raw
  )
  
  # Actual cell area, km2.
  area_current <- cellSize(
    current_raw,
    unit = "km"
  )
  
  
  for (j in seq_len(nrow(future_design))) {
    
    period_i <- future_design$period[j]
    scenario_i <- future_design$scenario[j]
    scenario_label_i <- future_design$scenario_label[j]
    
    future_file <- file.path(
      group_dir_i,
      future_design$future_name[j]
    )
    
    cat(
      "  ",
      period_i,
      " | ",
      scenario_label_i,
      "\n",
      sep = ""
    )
    
    future_raw <- rast(
      future_file
    )
    
    future_raw <- align_to_template(
      current_raw,
      future_raw,
      method = alignment_method
    )
    
    future_p_full <- to_prob01(
      future_raw,
      future_file
    )
    
    # Strict common domain:
    # current suitability + future suitability + DEM + latitude + area
    common_valid <- ifel(
      !is.na(current_p_full) &
        !is.na(future_p_full) &
        !is.na(dem_current) &
        !is.na(abs_lat_current) &
        !is.na(area_current),
      1,
      NA
    )
    
    current_p <- mask(
      current_p_full,
      common_valid
    )
    
    future_p <- mask(
      future_p_full,
      common_valid
    )
    
    dem_use <- mask(
      dem_current,
      common_valid
    )
    
    lat_use <- mask(
      abs_lat_current,
      common_valid
    )
    
    area_use <- mask(
      area_current,
      common_valid
    )
    
    
    # --------------------------------------------------------
    # Latitude metrics
    # --------------------------------------------------------
    
    lat_mean_current <- weighted_mean_raster(
      current_p,
      lat_use,
      area_use
    )
    
    lat_mean_future <- weighted_mean_raster(
      future_p,
      lat_use,
      area_use
    )
    
    lat_shift <- lat_mean_future -
      lat_mean_current
    
    lat_q_current <- weighted_quantile_raster(
      current_p,
      lat_use,
      area_use,
      probs = c(0.05, 0.50, 0.95)
    )
    
    lat_q_future <- weighted_quantile_raster(
      future_p,
      lat_use,
      area_use,
      probs = c(0.05, 0.50, 0.95)
    )
    
    lat_q05_shift <- lat_q_future["q5"] -
      lat_q_current["q5"]
    
    lat_q50_shift <- lat_q_future["q50"] -
      lat_q_current["q50"]
    
    lat_q95_shift <- lat_q_future["q95"] -
      lat_q_current["q95"]
    
    
    # --------------------------------------------------------
    # Elevation metrics
    # --------------------------------------------------------
    
    elev_mean_current <- weighted_mean_raster(
      current_p,
      dem_use,
      area_use
    )
    
    elev_mean_future <- weighted_mean_raster(
      future_p,
      dem_use,
      area_use
    )
    
    elev_shift <- elev_mean_future -
      elev_mean_current
    
    elev_q_current <- weighted_quantile_raster(
      current_p,
      dem_use,
      area_use,
      probs = c(0.05, 0.50, 0.95)
    )
    
    elev_q_future <- weighted_quantile_raster(
      future_p,
      dem_use,
      area_use,
      probs = c(0.05, 0.50, 0.95)
    )
    
    elev_q05_shift <- elev_q_future["q5"] -
      elev_q_current["q5"]
    
    elev_q50_shift <- elev_q_future["q50"] -
      elev_q_current["q50"]
    
    elev_q95_shift <- elev_q_future["q95"] -
      elev_q_current["q95"]
    
    
    # High-suitability core mean elevation.
    elev_core_current <- core_mean_elevation(
      current_p,
      dem_use,
      threshold = core_threshold
    )
    
    elev_core_future <- core_mean_elevation(
      future_p,
      dem_use,
      threshold = core_threshold
    )
    
    elev_core_shift <- elev_core_future -
      elev_core_current
    
    
    # --------------------------------------------------------
    # Save statistics
    # --------------------------------------------------------
    
    result_list[[result_id]] <- tibble(
      output_index = output_index_i,
      pole = pole_i,
      group_type = group_type_i,
      group_name = group_name_i,
      period = period_i,
      scenario = scenario_i,
      scenario_label = scenario_label_i,
      
      lat_mean_current = lat_mean_current,
      lat_mean_future = lat_mean_future,
      lat_shift = lat_shift,
      
      lat_q05_current = as.numeric(lat_q_current["q5"]),
      lat_q50_current = as.numeric(lat_q_current["q50"]),
      lat_q95_current = as.numeric(lat_q_current["q95"]),
      lat_q05_future = as.numeric(lat_q_future["q5"]),
      lat_q50_future = as.numeric(lat_q_future["q50"]),
      lat_q95_future = as.numeric(lat_q_future["q95"]),
      lat_q05_shift = as.numeric(lat_q05_shift),
      lat_q50_shift = as.numeric(lat_q50_shift),
      lat_q95_shift = as.numeric(lat_q95_shift),
      
      elev_mean_current = elev_mean_current,
      elev_mean_future = elev_mean_future,
      elev_shift = elev_shift,
      
      elev_q05_current = as.numeric(elev_q_current["q5"]),
      elev_q50_current = as.numeric(elev_q_current["q50"]),
      elev_q95_current = as.numeric(elev_q_current["q95"]),
      elev_q05_future = as.numeric(elev_q_future["q5"]),
      elev_q50_future = as.numeric(elev_q_future["q50"]),
      elev_q95_future = as.numeric(elev_q_future["q95"]),
      elev_q05_shift = as.numeric(elev_q05_shift),
      elev_q50_shift = as.numeric(elev_q50_shift),
      elev_q95_shift = as.numeric(elev_q95_shift),
      
      elev_core_current = elev_core_current,
      elev_core_future = elev_core_future,
      elev_core_shift = elev_core_shift
    )
    
    result_id <- result_id + 1L
    
    
    # --------------------------------------------------------
    # 2090s SSP5-8.5:
    # weighted profiles + exact spatial contribution maps
    # --------------------------------------------------------
    
    if (
      period_i == target_period &&
      scenario_i == target_scenario
    ) {
      
      profile_list[[profile_id]] <- bind_rows(
        
        make_weighted_profile(
          prob_r = current_p,
          covar_r = lat_use,
          area_r = area_use,
          bin_width = latitude_bin_width_deg,
          variable_name = "Absolute latitude",
          state_name = "Current",
          pole = pole_i,
          group_type = group_type_i,
          group_name = group_name_i
        ),
        
        make_weighted_profile(
          prob_r = future_p,
          covar_r = lat_use,
          area_r = area_use,
          bin_width = latitude_bin_width_deg,
          variable_name = "Absolute latitude",
          state_name = "2090s SSP5-8.5",
          pole = pole_i,
          group_type = group_type_i,
          group_name = group_name_i
        ),
        
        make_weighted_profile(
          prob_r = current_p,
          covar_r = dem_use,
          area_r = area_use,
          bin_width = elevation_bin_width_m,
          variable_name = "Elevation",
          state_name = "Current",
          pole = pole_i,
          group_type = group_type_i,
          group_name = group_name_i
        ),
        
        make_weighted_profile(
          prob_r = future_p,
          covar_r = dem_use,
          area_r = area_use,
          bin_width = elevation_bin_width_m,
          variable_name = "Elevation",
          state_name = "2090s SSP5-8.5",
          pole = pole_i,
          group_type = group_type_i,
          group_name = group_name_i
        )
      )
      
      profile_id <- profile_id + 1L
      
      
      # Exact spatial contribution to latitude shift.
      lat_contribution <- make_shift_contribution_raster(
        current_p = current_p,
        future_p = future_p,
        covar_r = lat_use,
        area_r = area_use
      )
      
      # Exact spatial contribution to elevation shift.
      elev_contribution <- make_shift_contribution_raster(
        current_p = current_p,
        future_p = future_p,
        covar_r = dem_use,
        area_r = area_use
      )
      
      group_contribution_dir <- file.path(
        contribution_raster_dir,
        pole_i,
        group_type_i,
        group_name_i
      )
      
      lat_contribution_file <- write_float_raster(
        lat_contribution,
        file.path(
          group_contribution_dir,
          "2090_SSP585_latitude_shift_contribution.tif"
        )
      )
      
      elev_contribution_file <- write_float_raster(
        elev_contribution,
        file.path(
          group_contribution_dir,
          "2090_SSP585_elevation_shift_contribution.tif"
        )
      )
      
      
      # Numerical verification:
      # sum of contribution raster must equal weighted-mean shift.
      lat_contribution_sum <- g_sum(
        lat_contribution
      )
      
      elev_contribution_sum <- g_sum(
        elev_contribution
      )
      
      cat(
        "    Latitude shift = ",
        round(lat_shift, 6),
        " deg | contribution sum = ",
        round(lat_contribution_sum, 6),
        "\n",
        sep = ""
      )
      
      cat(
        "    Elevation shift = ",
        round(elev_shift, 3),
        " m | contribution sum = ",
        round(elev_contribution_sum, 3),
        "\n",
        sep = ""
      )
      
      
      contribution_manifest_list[[manifest_id]] <- tibble(
        output_index = output_index_i,
        pole = pole_i,
        group_type = group_type_i,
        group_name = group_name_i,
        lat_contribution_file = lat_contribution_file,
        elev_contribution_file = elev_contribution_file,
        lat_shift = lat_shift,
        elev_shift = elev_shift
      )
      
      manifest_id <- manifest_id + 1L
    }
  }
}


# ============================================================
# 10. Export metric tables
# ============================================================

shift_summary_df <- bind_rows(
  result_list
)

profile_df <- bind_rows(
  profile_list
)

contribution_manifest_df <- bind_rows(
  contribution_manifest_list
)

write_csv(
  shift_summary_df,
  file.path(
    table_dir,
    "Fig6_all_latitude_elevation_shift_metrics.csv"
  ),
  na = ""
)

write_csv(
  profile_df,
  file.path(
    table_dir,
    "Fig6_2090s_SSP585_weighted_profiles.csv"
  ),
  na = ""
)

write_csv(
  contribution_manifest_df,
  file.path(
    table_dir,
    "Fig6_2090s_SSP585_shift_contribution_manifest.csv"
  ),
  na = ""
)


# ============================================================
# 11. Figure 6 summary:
#     latitude and elevation shifts across all future scenarios
#
# Changes in this version:
# - four SSP values for each group/period are explicitly connected by
#   grey dashed segments;
# - the scenario legend is a separate manual panel on the RIGHT,
#   so it is guaranteed to remain vertical and compact.
# ============================================================

plot_order <- groups %>%
  select(
    pole,
    group_type,
    group_name
  ) %>%
  mutate(
    group_label = paste(
      pole,
      group_name,
      sep = " | "
    ),
    group_label = factor(
      group_label,
      levels = rev(
        paste(
          pole,
          group_name,
          sep = " | "
        )
      )
    )
  )

scenario_levels <- c(
  "SSP1-2.6",
  "SSP2-4.5",
  "SSP3-7.0",
  "SSP5-8.5"
)

summary_plot_df <- shift_summary_df %>%
  left_join(
    plot_order,
    by = c(
      "pole",
      "group_type",
      "group_name"
    )
  ) %>%
  mutate(
    scenario_label = factor(
      scenario_label,
      levels = scenario_levels
    ),
    scenario_id = as.integer(
      scenario_label
    ),
    period = factor(
      period,
      levels = c(
        "2050s",
        "2090s"
      )
    )
  ) %>%
  arrange(
    period,
    group_label,
    scenario_id
  )

scenario_colours <- c(
  "SSP1-2.6" = "#4C78A8",
  "SSP2-4.5" = "#72B7B2",
  "SSP3-7.0" = "#F2CF5B",
  "SSP5-8.5" = "#E45756"
)


make_shift_panel <- function(
    data,
    metric,
    x_label,
    panel_title
) {
  
  segment_df <- data %>%
    mutate(
      metric_value = .data[[metric]]
    ) %>%
    arrange(
      period,
      group_label,
      scenario_id
    ) %>%
    group_by(
      period,
      group_label
    ) %>%
    mutate(
      metric_next = lead(
        metric_value
      )
    ) %>%
    ungroup() %>%
    filter(
      !is.na(metric_next)
    )
  
  ggplot(
    data,
    aes(
      x = .data[[metric]],
      y = group_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      colour = "grey65",
      linewidth = 0.30,
      linetype = "dashed"
    ) +
    
    # Explicit horizontal links between adjacent SSPs.
    geom_segment(
      data = segment_df,
      aes(
        x = metric_value,
        xend = metric_next,
        y = group_label,
        yend = group_label
      ),
      inherit.aes = FALSE,
      colour = "grey55",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    
    geom_point(
      aes(
        colour = scenario_label
      ),
      size = 1.8,
      show.legend = FALSE
    ) +
    
    facet_wrap(
      ~period,
      nrow = 1
    ) +
    
    scale_colour_manual(
      values = scenario_colours,
      limits = scenario_levels,
      drop = FALSE
    ) +
    
    labs(
      title = panel_title,
      x = x_label,
      y = NULL
    ) +
    
    theme_classic(
      base_size = base_font_pt,
      base_family = base_family
    ) +
    
    theme(
      plot.title = element_text(
        size = base_font_pt + 0.5,
        face = "bold",
        hjust = 0
      ),
      strip.text = element_text(
        size = base_font_pt,
        face = "bold"
      ),
      axis.text = element_text(
        size = base_font_pt,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = base_font_pt,
        face = "bold"
      ),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}


p_lat_shift <- make_shift_panel(
  data = summary_plot_df,
  metric = "lat_shift",
  x_label = expression(
    Delta * " absolute latitude (degrees)"
  ),
  panel_title = "Poleward redistribution"
)

p_elev_shift <- make_shift_panel(
  data = summary_plot_df,
  metric = "elev_shift",
  x_label = expression(
    Delta * " elevation (m)"
  ),
  panel_title = "Upslope redistribution"
)


# ------------------------------------------------------------
# Manual vertical scenario legend
# ------------------------------------------------------------

legend_df <- tibble(
  scenario_label = factor(
    scenario_levels,
    levels = scenario_levels
  ),
  x = 0,
  y = rev(
    seq_along(
      scenario_levels
    )
  )
)

p_scenario_legend <- ggplot(
  legend_df,
  aes(
    x = x,
    y = y
  )
) +
  annotate(
    "text",
    x = 0,
    y = 5.0,
    label = "Scenario",
    family = base_family,
    fontface = "bold",
    hjust = 0,
    size = base_font_pt / 2.35
  ) +
  geom_point(
    aes(
      colour = scenario_label
    ),
    size = 2.2,
    show.legend = FALSE
  ) +
  geom_text(
    aes(
      x = 0.30,
      label = scenario_label
    ),
    family = base_family,
    size = base_font_pt / 2.35,
    hjust = 0
  ) +
  scale_colour_manual(
    values = scenario_colours,
    limits = scenario_levels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      -0.15,
      1.85
    ),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(
      0.4,
      5.35
    ),
    expand = c(0, 0)
  ) +
  theme_void(
    base_size = base_font_pt,
    base_family = base_family
  )


p_summary_main <- patchwork::wrap_plots(
  p_lat_shift,
  p_elev_shift,
  ncol = 1,
  heights = c(
    1,
    1
  )
)

p_fig6_summary <- patchwork::wrap_plots(
  p_summary_main,
  p_scenario_legend,
  ncol = 2,
  widths = c(
    5.4,
    1.25
  )
) +
  patchwork::plot_annotation(
    title = "Climate-driven redistribution of climatically suitable habitats"
  )


summary_tiff <- file.path(
  summary_figure_dir,
  "Fig6_latitude_elevation_shifts.tiff"
)

summary_pdf <- file.path(
  summary_figure_dir,
  "Fig6_latitude_elevation_shifts.pdf"
)

summary_png <- file.path(
  summary_figure_dir,
  "Fig6_latitude_elevation_shifts.png"
)

ggsave(
  summary_tiff,
  p_fig6_summary,
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  summary_pdf,
  p_fig6_summary,
  device = grDevices::cairo_pdf,
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm"
)

ggsave(
  summary_png,
  p_fig6_summary,
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm",
  dpi = 600
)


# ============================================================
# 12. Weighted distribution profiles:
#     Current vs 2090s SSP5-8.5
#
# Each pole gets one 6-row × 2-column figure:
# left  = absolute-latitude profile
# right = elevation profile
#
# A semi-transparent band shows the absolute difference between
# current and future weighted distributions.
# ============================================================

profile_state_colours <- c(
  "Current" = "grey30",
  "2090s SSP5-8.5" = "#E66101"
)

profile_df <- profile_df %>%
  mutate(
    variable = factor(
      variable,
      levels = c(
        "Absolute latitude",
        "Elevation"
      )
    ),
    group_name = factor(
      group_name,
      levels = c(
        "Pale",
        "Bright",
        "Dark",
        "Crustose",
        "Foliose",
        "Fruticose"
      )
    )
  )


for (pole_i in c("Antarctic", "Arctic")) {
  
  profile_pole_df <- profile_df %>%
    filter(
      pole == pole_i
    )
  
  profile_band_df <- profile_pole_df %>%
    select(
      pole,
      group_name,
      variable,
      bin,
      state,
      weighted_fraction
    ) %>%
    tidyr::pivot_wider(
      names_from = state,
      values_from = weighted_fraction
    ) %>%
    mutate(
      Current = ifelse(
        is.na(Current),
        0,
        Current
      ),
      `2090s SSP5-8.5` = ifelse(
        is.na(`2090s SSP5-8.5`),
        0,
        `2090s SSP5-8.5`
      ),
      ymin = pmin(
        Current,
        `2090s SSP5-8.5`
      ),
      ymax = pmax(
        Current,
        `2090s SSP5-8.5`
      ),
      diff_value =
        `2090s SSP5-8.5` -
        Current
    )
  
  p_profile <- ggplot() +
    
    # Difference band.
    geom_ribbon(
      data = profile_band_df,
      aes(
        x = bin,
        ymin = ymin,
        ymax = ymax
      ),
      fill = "#E66101",
      alpha = 0.16
    ) +
    
    # Current/future profile lines.
    geom_line(
      data = profile_pole_df,
      aes(
        x = bin,
        y = weighted_fraction,
        colour = state
      ),
      linewidth = 0.55
    ) +
    
    facet_grid(
      group_name ~ variable,
      scales = "free_x"
    ) +
    
    scale_colour_manual(
      values = profile_state_colours,
      breaks = c(
        "2090s SSP5-8.5",
        "Current"
      )
    ) +
    
    labs(
      title = paste0(
        pole_i,
        " suitability-weighted redistribution under SSP5-8.5 in the 2090s"
      ),
      subtitle = "Shaded band indicates the difference between current and future weighted distributions.",
      x = NULL,
      y = "Fraction of WSA",
      colour = NULL
    ) +
    
    theme_classic(
      base_size = base_font_pt,
      base_family = base_family
    ) +
    
    theme(
      plot.title = element_text(
        size = base_font_pt + 1,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = base_font_pt,
        hjust = 0
      ),
      strip.text = element_text(
        size = base_font_pt,
        face = "bold"
      ),
      axis.text = element_text(
        size = base_font_pt,
        colour = "black"
      ),
      axis.title.y = element_text(
        size = base_font_pt,
        face = "bold"
      ),
      legend.position = "top",
      legend.text = element_text(
        size = base_font_pt
      )
    )
  
  profile_stem <- paste0(
    "Fig6_",
    pole_i,
    "_2090s_SSP585_weighted_latitude_elevation_profiles"
  )
  
  ggsave(
    file.path(
      profile_dir,
      paste0(
        profile_stem,
        ".tiff"
      )
    ),
    p_profile,
    width = profile_width_mm,
    height = profile_height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(
      profile_dir,
      paste0(
        profile_stem,
        ".pdf"
      )
    ),
    p_profile,
    device = grDevices::cairo_pdf,
    width = profile_width_mm,
    height = profile_height_mm,
    units = "mm"
  )
  
  ggsave(
    file.path(
      profile_dir,
      paste0(
        profile_stem,
        ".png"
      )
    ),
    p_profile,
    width = profile_width_mm,
    height = profile_height_mm,
    units = "mm",
    dpi = 600
  )
}


# ============================================================
# 13. Exact spatial contribution maps:
#     2090s SSP5-8.5 only
#
# Left:
# positive red = contribution to poleward shift
# negative blue = contribution to equatorward shift
#
# Right:
# positive red = contribution to upslope shift
# negative blue = contribution to downslope shift
# ============================================================

graphics.off()

for (i in seq_len(nrow(contribution_manifest_df))) {
  
  pole_i <- contribution_manifest_df$pole[i]
  group_type_i <- contribution_manifest_df$group_type[i]
  group_name_i <- contribution_manifest_df$group_name[i]
  output_index_i <- contribution_manifest_df$output_index[i]
  
  outline_df <- outline_df_list[[pole_i]]
  
  p_lat_contrib <- make_contribution_map(
    raster_file = contribution_manifest_df$lat_contribution_file[i],
    pole = pole_i,
    outline_df = outline_df,
    panel_title = "Contribution to latitude shift",
    legend_title = expression(Delta * "|" * latitude * "|" * " contribution")
  )
  
  p_elev_contrib <- make_contribution_map(
    raster_file = contribution_manifest_df$elev_contribution_file[i],
    pole = pole_i,
    outline_df = outline_df,
    panel_title = "Contribution to elevation shift",
    legend_title = expression(Delta * "elevation contribution")
  )
  
  combined_contribution <- patchwork::wrap_plots(
    p_lat_contrib,
    p_elev_contrib,
    widths = c(1, 1)
  ) +
    patchwork::plot_annotation(
      title = paste0(
        pole_i,
        " | ",
        group_type_i,
        " | ",
        group_name_i,
        " | 2090s SSP5-8.5"
      ),
      subtitle = paste0(
        "Latitude shift = ",
        sprintf(
          "%+.3f°",
          contribution_manifest_df$lat_shift[i]
        ),
        "; elevation shift = ",
        sprintf(
          "%+.1f m",
          contribution_manifest_df$elev_shift[i]
        ),
        ". Red areas contribute positively to poleward/upslope redistribution; ",
        "blue areas contribute negatively."
      )
    )
  
  out_stem <- paste(
    output_index_i,
    pole_i,
    gsub(
      " ",
      "_",
      group_type_i
    ),
    group_name_i,
    "2090s_SSP585_lat_elev_contribution",
    sep = "_"
  )
  
  ggsave(
    file.path(
      contribution_figure_dir,
      paste0(out_stem, ".tiff")
    ),
    combined_contribution,
    width = contribution_width_mm,
    height = contribution_height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    file.path(
      contribution_figure_dir,
      paste0(out_stem, ".pdf")
    ),
    combined_contribution,
    device = grDevices::cairo_pdf,
    width = contribution_width_mm,
    height = contribution_height_mm,
    units = "mm"
  )
  
  ggsave(
    file.path(
      contribution_figure_dir,
      paste0(out_stem, ".png")
    ),
    combined_contribution,
    width = contribution_width_mm,
    height = contribution_height_mm,
    units = "mm",
    dpi = 600
  )
  
  graphics.off()
}


# ============================================================
# 14. Key 2090s SSP5-8.5 table for manuscript checking
# ============================================================

key_2090_ssp585 <- shift_summary_df %>%
  filter(
    period == target_period,
    scenario == target_scenario
  ) %>%
  select(
    pole,
    group_type,
    group_name,
    
    lat_mean_current,
    lat_mean_future,
    lat_shift,
    lat_q05_shift,
    lat_q50_shift,
    lat_q95_shift,
    
    elev_mean_current,
    elev_mean_future,
    elev_shift,
    elev_q05_shift,
    elev_q50_shift,
    elev_q95_shift,
    
    elev_core_current,
    elev_core_future,
    elev_core_shift
  )

write_csv(
  key_2090_ssp585,
  file.path(
    table_dir,
    "Fig6_2090s_SSP585_key_metrics.csv"
  ),
  na = ""
)

print(
  key_2090_ssp585
)


# ============================================================
# 15. Completion
# ============================================================

cat("\n==================================================\n")
cat("Completed latitude/elevation redistribution analysis.\n")
cat("All scenarios table:\n  ",
    file.path(table_dir, "Fig6_all_latitude_elevation_shift_metrics.csv"),
    "\n", sep = "")
cat("Key 2090s SSP5-8.5 table:\n  ",
    file.path(table_dir, "Fig6_2090s_SSP585_key_metrics.csv"),
    "\n", sep = "")
cat("Summary figure:\n  ", summary_png, "\n", sep = "")
cat("Weighted profiles:\n  ", profile_dir, "\n", sep = "")
cat("Spatial contribution maps:\n  ", contribution_figure_dir, "\n", sep = "")
cat("Contribution rasters:\n  ", contribution_raster_dir, "\n", sep = "")
cat("==================================================\n")















# ============================================================
# ============================================================
# ============================================================
# Fig. 6: latitude and elevation redistribution
# 1) Summary figure: delta absolute latitude + delta elevation
# 2) Optional trajectory figure: multidimensional redistribution
# 3) Antarctic weighted profile figure (Current vs 2090s SSP5-8.5)
# 4) Arctic weighted profile figure (Current vs 2090s SSP5-8.5)
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "forcats",
  "patchwork",
  "ggrepel",
  "scales",
  "stringr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
lapply(required_packages, library, character.only = TRUE)

# -----------------------------
# 1. Paths
# -----------------------------
metrics_file  <- "F:/1 MY PROJECTS/Climate Change/分析过程/_Fig6_latitude_elevation_redistribution/tables/Fig6_all_latitude_elevation_shift_metrics.csv"
profiles_file <- "F:/1 MY PROJECTS/Climate Change/分析过程/_Fig6_latitude_elevation_redistribution/tables/Fig6_2090s_SSP585_weighted_profiles.csv"

out_dir <- "F:/1 MY PROJECTS/Climate Change/分析过程/_Fig6_latitude_elevation_redistribution/Fig6_outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 2. Global plotting settings
# -----------------------------
base_family   <- "Arial"
base_font_pt  <- 9

scenario_levels <- c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
scenario_cols <- c(
  "SSP1-2.6" = "#4C78A8",
  "SSP2-4.5" = "#72B7B2",
  "SSP3-7.0" = "#EAC54F",
  "SSP5-8.5" = "#E45756"
)

group_order <- c("Pale", "Bright", "Dark", "Crustose", "Foliose", "Fruticose")
pole_order  <- c("Antarctic", "Arctic")
period_order <- c("2050s", "2090s")

# y 位置：南极 6 行 + 1 行空白 + 北极 6 行
group_layout <- tibble::tibble(
  pole       = c(rep("Antarctic", 6), rep("Arctic", 6)),
  group_name = rep(group_order, 2),
  y          = c(12, 11, 10, 9, 8, 7, 5, 4, 3, 2, 1, 0),
  group_label = paste(pole, group_name, sep = " | ")
)

# 通用主题
theme_fig6 <- function() {
  theme_bw(base_size = base_font_pt, base_family = base_family) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = base_font_pt, face = "bold"),
      axis.text  = element_text(size = base_font_pt),
      strip.text = element_text(size = base_font_pt, face = "bold"),
      plot.title = element_text(size = base_font_pt + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_font_pt, hjust = 0),
      legend.title = element_text(size = base_font_pt, face = "bold"),
      legend.text  = element_text(size = base_font_pt),
      legend.position = "right",
      legend.direction = "vertical"
    )
}

# 保存函数
save_plot_multi <- function(plot_obj, out_stem, width_mm, height_mm) {
  ggsave(
    filename = file.path(out_dir, paste0(out_stem, ".png")),
    plot = plot_obj,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600
  )
  ggsave(
    filename = file.path(out_dir, paste0(out_stem, ".pdf")),
    plot = plot_obj,
    width = width_mm,
    height = height_mm,
    units = "mm",
    device = cairo_pdf
  )
}

# -----------------------------
# 3. Read and prepare data
# -----------------------------
metrics_df <- readr::read_csv(metrics_file, show_col_types = FALSE)
profiles_df <- readr::read_csv(profiles_file, show_col_types = FALSE)

metrics_df <- metrics_df %>%
  mutate(
    pole = factor(pole, levels = pole_order),
    group_name = factor(group_name, levels = group_order),
    period = factor(period, levels = period_order),
    scenario_label = factor(scenario_label, levels = scenario_levels),
    group_label = paste(as.character(pole), as.character(group_name), sep = " | ")
  ) %>%
  left_join(group_layout, by = c("group_label", "group_name", "pole"))

profiles_df <- profiles_df %>%
  mutate(
    pole = factor(pole, levels = pole_order),
    group_name = factor(group_name, levels = group_order),
    variable = factor(variable, levels = c("Absolute latitude", "Elevation")),
    state = factor(state, levels = c("Current", "2090s SSP5-8.5"))
  )

# -----------------------------
# 4. Fig. 6 summary figure
#    a) Delta absolute latitude
#    b) Delta elevation
# -----------------------------

# connector segments：同一组、同一时期，4 个 SSP 之间用灰色虚线横向连接
lat_segments <- metrics_df %>%
  group_by(period, pole, group_name, y) %>%
  summarise(
    xmin = min(lat_shift, na.rm = TRUE),
    xmax = max(lat_shift, na.rm = TRUE),
    .groups = "drop"
  )

elev_segments <- metrics_df %>%
  group_by(period, pole, group_name, y) %>%
  summarise(
    xmin = min(elev_shift, na.rm = TRUE),
    xmax = max(elev_shift, na.rm = TRUE),
    .groups = "drop"
  )

# 面板 a：纬度变化
p_lat <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5, linetype = "22") +
  geom_segment(
    data = lat_segments,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 0.5,
    linetype = "22"
  ) +
  geom_point(
    data = metrics_df,
    aes(x = lat_shift, y = y, colour = scenario_label),
    size = 2.8
  ) +
  geom_hline(yintercept = 6, colour = "grey80", linewidth = 0.5) +
  facet_wrap(~ period, nrow = 1) +
  scale_y_continuous(
    breaks = group_layout$y,
    labels = group_layout$group_label,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_colour_manual(values = scenario_cols, drop = FALSE) +
  labs(
    title = "a  Poleward redistribution",
    x = expression(Delta * " absolute latitude (degrees)"),
    y = NULL,
    colour = "Scenario"
  ) +
  coord_cartesian(clip = "off") +
  theme_fig6() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.text.y = element_text(hjust = 1)
  )

# 面板 b：海拔变化
p_elev <- ggplot() +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.5, linetype = "22") +
  geom_segment(
    data = elev_segments,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 0.5,
    linetype = "22"
  ) +
  geom_point(
    data = metrics_df,
    aes(x = elev_shift, y = y, colour = scenario_label),
    size = 2.8
  ) +
  geom_hline(yintercept = 6, colour = "grey80", linewidth = 0.5) +
  facet_wrap(~ period, nrow = 1) +
  scale_y_continuous(
    breaks = group_layout$y,
    labels = group_layout$group_label,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_colour_manual(values = scenario_cols, drop = FALSE) +
  labs(
    title = "b  Upslope redistribution",
    x = expression(Delta * " elevation (m)"),
    y = NULL,
    colour = "Scenario"
  ) +
  coord_cartesian(clip = "off") +
  theme_fig6() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "white", colour = "black"),
    axis.text.y = element_text(hjust = 1)
  )

# ------------------------------------------------------------
# Put the legend on the right before combining the two panels
# Compatible with older patchwork versions
# ------------------------------------------------------------

p_lat <- p_lat +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  )

p_elev <- p_elev +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  )


fig6_summary <- (p_lat / p_elev) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(1, 1)
  ) +
  patchwork::plot_annotation(
    title = "Climate-driven redistribution of climatically suitable habitats",
    subtitle = paste0(
      "Positive \u0394 absolute latitude indicates poleward redistribution; ",
      "positive \u0394 elevation indicates upslope redistribution."
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family,
        size = base_font_pt + 3,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        family = base_family,
        size = base_font_pt,
        hjust = 0
      )
    )
  )

save_plot_multi(fig6_summary, "Fig6_summary_latitude_elevation_shift", 180, 240)

# -----------------------------
# 5. Optional Fig. 6c
#    Multidimensional redistribution trajectory
# -----------------------------
# 每个 group 在每个 pole × period 内，连接 4 个 SSP
traj_df <- metrics_df %>%
  arrange(pole, period, group_name, scenario_label)

traj_labels <- traj_df %>%
  filter(scenario_label == "SSP5-8.5")

p_traj <- ggplot(traj_df, aes(x = lat_shift, y = elev_shift, colour = scenario_label)) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = "22", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = "22", linewidth = 0.5) +
  geom_path(
    aes(group = interaction(pole, period, group_name)),
    colour = "grey55",
    linewidth = 0.5,
    linetype = "22",
    show.legend = FALSE
  ) +
  geom_point(size = 2.7) +
  ggrepel::geom_text_repel(
    data = traj_labels,
    aes(label = group_name),
    size = 2.6,
    family = base_family,
    colour = "black",
    min.segment.length = 0,
    box.padding = 0.18,
    point.padding = 0.1,
    show.legend = FALSE
  ) +
  facet_grid(pole ~ period) +
  scale_colour_manual(values = scenario_cols, drop = FALSE) +
  labs(
    title = "Multidimensional redistribution trajectory",
    x = expression(Delta * " absolute latitude (degrees)"),
    y = expression(Delta * " elevation (m)"),
    colour = "Scenario"
  ) +
  theme_fig6() +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    legend.position = "right"
  )

save_plot_multi(p_traj, "Fig6c_multidimensional_redistribution_trajectory", 180, 140)

# -----------------------------
# 6. Weighted profiles
#    Current vs 2090s SSP5-8.5
#    Add transparent ribbons for differences
# -----------------------------

profiles_wide <- profiles_df %>%
  select(bin, weighted_fraction, variable, state, pole, group_name) %>%
  tidyr::pivot_wider(
    names_from = state,
    values_from = weighted_fraction
  ) %>%
  mutate(
    diff = `2090s SSP5-8.5` - Current,
    diff_class = ifelse(diff >= 0, "Increase", "Decrease"),
    ymin = pmin(Current, `2090s SSP5-8.5`, na.rm = TRUE),
    ymax = pmax(Current, `2090s SSP5-8.5`, na.rm = TRUE)
  )

make_profile_plot <- function(pole_use) {
  
  line_df <- profiles_df %>%
    filter(pole == pole_use) %>%
    mutate(
      state = factor(state, levels = c("2090s SSP5-8.5", "Current"))
    )
  
  ribbon_df <- profiles_wide %>%
    filter(pole == pole_use)
  
  ggplot() +
    geom_ribbon(
      data = ribbon_df,
      aes(x = bin, ymin = ymin, ymax = ymax, fill = diff_class),
      alpha = 0.20,
      colour = NA
    ) +
    geom_line(
      data = line_df,
      aes(x = bin, y = weighted_fraction, colour = state),
      linewidth = 0.85
    ) +
    facet_grid(
      rows = vars(group_name),
      cols = vars(variable),
      scales = "free_x",
      switch = "y"
    ) +
    scale_colour_manual(
      values = c(
        "2090s SSP5-8.5" = "#E66101",
        "Current" = "grey35"
      ),
      breaks = c("2090s SSP5-8.5", "Current"),
      name = NULL
    ) +
    scale_fill_manual(
      values = c(
        "Increase" = "#3C78D8",
        "Decrease" = "#E66101"
      ),
      name = "Difference"
    ) +
    labs(
      title = paste0(pole_use, " suitability-weighted redistribution under SSP5-8.5 in the 2090s"),
      x = NULL,
      y = "Fraction of WSA"
    ) +
    theme_fig6() +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text.y.right = element_text(angle = 270, face = "bold"),
      legend.position = "top",
      legend.direction = "horizontal",
      panel.spacing = unit(0.6, "lines")
    )
}

p_profile_ant <- make_profile_plot("Antarctic")
p_profile_arc <- make_profile_plot("Arctic")

save_plot_multi(p_profile_ant, "Fig6_Antarctic_2090s_SSP585_weighted_profiles", 180, 250)
save_plot_multi(p_profile_arc, "Fig6_Arctic_2090s_SSP585_weighted_profiles", 180, 250)

# -----------------------------
# 7. Export source tables used in plotting (optional)
# -----------------------------
readr::write_csv(metrics_df, file.path(out_dir, "Fig6_metrics_used_in_plotting.csv"))
readr::write_csv(profiles_wide, file.path(out_dir, "Fig6_profiles_wide_for_ribbon.csv"))

cat("\nAll Fig. 6 outputs have been saved to:\n", out_dir, "\n")












# ============================================================
# Multidimensional redistribution trajectory
# Clean version
#
# x = lat_shift
# y = elev_shift
#
# Colour = SSP scenario
# Shape  = period
# Line   = same group + same SSP, connecting 2050s -> 2090s
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)


# ============================================================
# Fig. 6c
# Multidimensional redistribution of climatically suitable habitats
#
# Layout:
#   Antarctic | Color type    Antarctic | Growth form
#   Arctic    | Color type    Arctic    | Growth form
#
# Encoding:
#   Colour = SSP scenario
#   Circle = 2050s
#   Triangle = 2090s
#   Line = same functional group + same SSP, 2050s -> 2090s
#
# Arctic y-axis:      -50 to 50 m
# Antarctic y-axis:   -50 to 220 m
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(grid)


# ============================================================
# 1. Paths
# ============================================================

metrics_file <- paste0(
  "F:/1 MY PROJECTS/Climate Change/分析过程/",
  "_Fig6_latitude_elevation_redistribution/",
  "tables/Fig6_all_latitude_elevation_shift_metrics.csv"
)

trajectory_figure_dir <- paste0(
  "F:/1 MY PROJECTS/Climate Change/分析过程/",
  "_Fig6_latitude_elevation_redistribution/",
  "trajectory_figure"
)

dir.create(
  trajectory_figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Figure settings
# ============================================================

base_family <- "Arial"
base_font_pt <- 7

figure_width_mm  <- 180
figure_height_mm <- 145


scenario_levels <- c(
  "SSP1-2.6",
  "SSP2-4.5",
  "SSP3-7.0",
  "SSP5-8.5"
)


scenario_cols <- c(
  "SSP1-2.6" = "#4C78A8",
  "SSP2-4.5" = "#72B7B2",
  "SSP3-7.0" = "#E3BA44",
  "SSP5-8.5" = "#E45756"
)


period_shapes <- c(
  "2050s" = 16,   # circle
  "2090s" = 17    # triangle
)


group_levels <- c(
  "Pale",
  "Bright",
  "Dark",
  "Crustose",
  "Foliose",
  "Fruticose"
)


# ============================================================
# 3. Re-read original metrics table
# ============================================================

metrics_df <- readr::read_csv(
  metrics_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# Check required columns
# ------------------------------------------------------------

required_columns <- c(
  "pole",
  "group_type",
  "group_name",
  "period",
  "scenario_label",
  "lat_shift",
  "elev_shift"
)

missing_columns <- setdiff(
  required_columns,
  names(metrics_df)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ============================================================
# 4. Standardize data
# ============================================================

traj_df <- metrics_df %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Preserve original period for checking
    # --------------------------------------------------------
    
    period_raw = trimws(
      as.character(period)
    ),
    
    
    # --------------------------------------------------------
    # Standardize period
    # --------------------------------------------------------
    
    period = case_when(
      
      grepl(
        "^2050",
        period_raw
      ) ~ "2050s",
      
      grepl(
        "^2090",
        period_raw
      ) ~ "2090s",
      
      TRUE ~ NA_character_
    ),
    
    
    # --------------------------------------------------------
    # Standardize pole
    # --------------------------------------------------------
    
    pole = case_when(
      
      grepl(
        "^Antar",
        as.character(pole),
        ignore.case = TRUE
      ) ~ "Antarctic",
      
      grepl(
        "^Arctic",
        as.character(pole),
        ignore.case = TRUE
      ) ~ "Arctic",
      
      TRUE ~ as.character(pole)
    ),
    
    
    # --------------------------------------------------------
    # Reconstruct functional dimension from group name
    # This avoids any inconsistencies in group_type
    # --------------------------------------------------------
    
    group_type = case_when(
      
      as.character(group_name) %in% c(
        "Pale",
        "Bright",
        "Dark"
      ) ~ "Color type",
      
      as.character(group_name) %in% c(
        "Crustose",
        "Foliose",
        "Fruticose"
      ) ~ "Growth form",
      
      TRUE ~ as.character(group_type)
    ),
    
    
    # --------------------------------------------------------
    # Factors
    # --------------------------------------------------------
    
    pole = factor(
      pole,
      levels = c(
        "Antarctic",
        "Arctic"
      )
    ),
    
    group_type = factor(
      group_type,
      levels = c(
        "Color type",
        "Growth form"
      )
    ),
    
    group_name = factor(
      as.character(group_name),
      levels = group_levels
    ),
    
    period = factor(
      period,
      levels = c(
        "2050s",
        "2090s"
      )
    ),
    
    scenario_label = factor(
      as.character(scenario_label),
      levels = scenario_levels
    )
  ) %>%
  
  filter(
    !is.na(pole),
    !is.na(group_type),
    !is.na(group_name),
    !is.na(period),
    !is.na(scenario_label),
    is.finite(lat_shift),
    is.finite(elev_shift)
  ) %>%
  
  arrange(
    pole,
    group_type,
    group_name,
    scenario_label,
    period
  )


# ============================================================
# 5. Data checks
# ============================================================

cat("\n========================================\n")
cat("Trajectory data check\n")
cat("========================================\n")


print(
  traj_df %>%
    count(
      pole,
      group_type,
      period,
      scenario_label
    )
)


# Every pole × group type × period × SSP should contain 3 groups.

trajectory_count_check <- traj_df %>%
  count(
    pole,
    group_type,
    group_name,
    scenario_label
  )


if (any(trajectory_count_check$n != 2)) {
  
  print(
    trajectory_count_check %>%
      filter(
        n != 2
      )
  )
  
  stop(
    "Some functional-group/scenario trajectories ",
    "do not contain exactly two periods."
  )
}


cat(
  "\nAll trajectories contain exactly ",
  "two periods: 2050s and 2090s.\n"
)


# ============================================================
# 6. Build explicit 2050s -> 2090s line segments
# ============================================================
#
# This is more robust than geom_path().
# Each row below represents ONE trajectory:
#
# group + SSP:
#     2050s circle -------- 2090s triangle
# ============================================================

segment_df <- traj_df %>%
  
  select(
    pole,
    group_type,
    group_name,
    scenario_label,
    period,
    lat_shift,
    elev_shift
  ) %>%
  
  pivot_wider(
    names_from = period,
    values_from = c(
      lat_shift,
      elev_shift
    ),
    names_sep = "_"
  ) %>%
  
  filter(
    is.finite(lat_shift_2050s),
    is.finite(lat_shift_2090s),
    is.finite(elev_shift_2050s),
    is.finite(elev_shift_2090s)
  )


cat(
  "Number of 2050s -> 2090s trajectories: ",
  nrow(segment_df),
  "\n",
  sep = ""
)


# Expected:
# 2 poles × 6 groups × 4 SSP = 48 trajectories

if (nrow(segment_df) != 48) {
  
  warning(
    "Expected 48 trajectories, but found ",
    nrow(segment_df),
    "."
  )
}

# ============================================================
# 6.1 Build scenario-to-scenario segments
#
# Within the SAME:
#   pole
#   functional group
#   period
#
# connect:
# SSP1-2.6 -> SSP2-4.5 -> SSP3-7.0 -> SSP5-8.5
#
# These will be drawn as light-grey solid lines.
# ============================================================

scenario_segment_df <- traj_df %>%
  
  mutate(
    scenario_id = as.integer(
      scenario_label
    )
  ) %>%
  
  arrange(
    pole,
    group_type,
    group_name,
    period,
    scenario_id
  ) %>%
  
  group_by(
    pole,
    group_type,
    group_name,
    period
  ) %>%
  
  mutate(
    lat_shift_next = lead(
      lat_shift
    ),
    
    elev_shift_next = lead(
      elev_shift
    )
  ) %>%
  
  ungroup() %>%
  
  filter(
    !is.na(lat_shift_next),
    !is.na(elev_shift_next)
  )


# ============================================================
# 7. Labels
# ============================================================
#
# Only label:
# 2090s | SSP5-8.5
#
# Each panel therefore contains only three group labels.
# ============================================================

label_df <- traj_df %>%
  
  filter(
    period == "2090s",
    scenario_label == "SSP5-8.5"
  )


if (nrow(label_df) != 12) {
  
  warning(
    "Expected 12 endpoint labels, but found ",
    nrow(label_df),
    "."
  )
}


# ============================================================
# 8. Calculate common x-axis limits within each pole
# ============================================================
#
# Color type and Growth form share the same x-range
# within Antarctic or Arctic, improving comparability.
# ============================================================

make_x_limits <- function(
    data,
    padding_fraction = 0.06
) {
  
  xmin <- min(
    c(
      0,
      data$lat_shift
    ),
    na.rm = TRUE
  )
  
  xmax <- max(
    c(
      0,
      data$lat_shift
    ),
    na.rm = TRUE
  )
  
  xrange <- xmax - xmin
  
  if (
    !is.finite(xrange) ||
    xrange <= 0
  ) {
    xrange <- 1
  }
  
  c(
    xmin - padding_fraction * xrange,
    xmax + padding_fraction * xrange
  )
}


ant_x_limits <- make_x_limits(
  traj_df %>%
    filter(
      pole == "Antarctic"
    )
)


arc_x_limits <- make_x_limits(
  traj_df %>%
    filter(
      pole == "Arctic"
    )
)


cat(
  "\nAntarctic x limits: ",
  paste(
    round(
      ant_x_limits,
      3
    ),
    collapse = " to "
  ),
  "\n",
  sep = ""
)


cat(
  "Arctic x limits: ",
  paste(
    round(
      arc_x_limits,
      3
    ),
    collapse = " to "
  ),
  "\n",
  sep = ""
)


# ============================================================
# 9. Function for one panel
# ============================================================

make_trajectory_panel <- function(
    pole_use,
    group_type_use,
    panel_label,
    x_limits,
    y_limits,
    show_legend = FALSE,
    show_x_title = TRUE,
    show_y_title = TRUE
) {
  
  # ----------------------------------------------------------
  # Point data
  # ----------------------------------------------------------
  
  dat <- traj_df %>%
    filter(
      pole == pole_use,
      group_type == group_type_use
    )
  
  
  # ----------------------------------------------------------
  # Temporal segments:
  # same group + same SSP
  # 2050s -> 2090s
  # ----------------------------------------------------------
  
  seg_dat <- segment_df %>%
    filter(
      pole == pole_use,
      group_type == group_type_use
    )
  
  
  # ----------------------------------------------------------
  # Scenario segments:
  # same group + same period
  # SSP1-2.6 -> SSP2-4.5 -> SSP3-7.0 -> SSP5-8.5
  # ----------------------------------------------------------
  
  scenario_seg_dat <- scenario_segment_df %>%
    filter(
      pole == pole_use,
      group_type == group_type_use
    )
  
  
  # ----------------------------------------------------------
  # Labels:
  # 2090s SSP5-8.5 only
  # ----------------------------------------------------------
  
  lab_dat <- label_df %>%
    filter(
      pole == pole_use,
      group_type == group_type_use
    )
  
  
  # ==========================================================
  # Plot
  # ==========================================================
  
  p <- ggplot() +
    
    # --------------------------------------------------------
  # Zero reference lines
  # --------------------------------------------------------
  
  geom_vline(
    xintercept = 0,
    colour = "grey72",
    linewidth = 0.30,
    linetype = "dashed"
  ) +
    
    geom_hline(
      yintercept = 0,
      colour = "grey72",
      linewidth = 0.30,
      linetype = "dashed"
    ) +
    
    
    # --------------------------------------------------------
  # Scenario-intensity connections
  #
  # Same functional group + same period:
  #
  # SSP1-2.6 -> SSP2-4.5 ->
  # SSP3-7.0 -> SSP5-8.5
  #
  # Light-grey solid lines
  # --------------------------------------------------------
  
  geom_segment(
    data = scenario_seg_dat,
    
    aes(
      x = lat_shift,
      y = elev_shift,
      
      xend = lat_shift_next,
      yend = elev_shift_next
    ),
    
    inherit.aes = FALSE,
    
    colour = "grey70",
    linewidth = 0.38,
    linetype = "solid",
    alpha = 0.90,
    lineend = "round"
  ) +
    
    
    # --------------------------------------------------------
  # Temporal connections
  #
  # Same functional group + same SSP:
  # 2050s -> 2090s
  #
  # Colour = scenario
  # --------------------------------------------------------
  
  geom_segment(
    data = seg_dat,
    
    aes(
      x = lat_shift_2050s,
      y = elev_shift_2050s,
      
      xend = lat_shift_2090s,
      yend = elev_shift_2090s,
      
      colour = scenario_label
    ),
    
    linewidth = 0.55,
    alpha = 0.78,
    lineend = "round"
  ) +
    
    
    # --------------------------------------------------------
  # Points
  #
  # colour = SSP
  # shape  = period
  # --------------------------------------------------------
  
  geom_point(
    data = dat,
    
    aes(
      x = lat_shift,
      y = elev_shift,
      colour = scenario_label,
      shape = period
    ),
    
    size = 2.0,
    stroke = 0.20
  ) +
    
    
    # --------------------------------------------------------
  # Labels:
  # SSP5-8.5 in the 2090s
  # --------------------------------------------------------
  
  ggrepel::geom_text_repel(
    data = lab_dat,
    
    aes(
      x = lat_shift,
      y = elev_shift,
      label = group_name
    ),
    
    inherit.aes = FALSE,
    
    family = base_family,
    size = 2.35,
    
    colour = "black",
    
    segment.colour = "grey55",
    segment.linewidth = 0.22,
    
    box.padding = 0.18,
    point.padding = 0.12,
    
    min.segment.length = 0,
    max.overlaps = Inf,
    
    seed = 123,
    
    show.legend = FALSE
  ) +
    
    
    # --------------------------------------------------------
  # Scenario colours
  # --------------------------------------------------------
  
  scale_colour_manual(
    values = scenario_cols,
    
    breaks = scenario_levels,
    limits = scenario_levels,
    
    drop = FALSE,
    
    name = "Scenario"
  ) +
    
    
    # --------------------------------------------------------
  # Period shapes
  # --------------------------------------------------------
  
  scale_shape_manual(
    values = period_shapes,
    
    breaks = c(
      "2050s",
      "2090s"
    ),
    
    limits = c(
      "2050s",
      "2090s"
    ),
    
    drop = FALSE,
    
    name = "Period"
  ) +
    
    
    # --------------------------------------------------------
  # Coordinates
  # --------------------------------------------------------
  
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits,
    clip = "off"
  ) +
    
    
    # --------------------------------------------------------
  # Axis labels
  # --------------------------------------------------------
  
  labs(
    title = panel_label,
    
    x = if (show_x_title) {
      expression(
        Delta * " absolute latitude (degrees)"
      )
    } else {
      NULL
    },
    
    y = if (show_y_title) {
      expression(
        Delta * " elevation (m)"
      )
    } else {
      NULL
    }
  ) +
    
    
    # --------------------------------------------------------
  # Theme
  # --------------------------------------------------------
  
  theme_classic(
    base_size = base_font_pt,
    base_family = base_family
  ) +
    
    theme(
      
      plot.title = element_text(
        family = base_family,
        size = base_font_pt + 0.5,
        face = "bold",
        hjust = 0
      ),
      
      axis.title = element_text(
        family = base_family,
        size = base_font_pt,
        face = "bold"
      ),
      
      axis.text = element_text(
        family = base_family,
        size = base_font_pt,
        colour = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.30,
        colour = "black"
      ),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.30
      ),
      
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      },
      
      legend.direction = "vertical",
      legend.box = "vertical",
      
      legend.title = element_text(
        family = base_family,
        size = base_font_pt,
        face = "bold"
      ),
      
      legend.text = element_text(
        family = base_family,
        size = base_font_pt
      ),
      
      plot.margin = margin(
        t = 3,
        r = 9,
        b = 3,
        l = 3,
        unit = "pt"
      )
    )
  
  
  # ----------------------------------------------------------
  # Legend
  # ----------------------------------------------------------
  
  if (show_legend) {
    
    p <- p +
      
      guides(
        
        colour = guide_legend(
          order = 1,
          title = "Scenario",
          override.aes = list(
            shape = 16,
            size = 2
          )
        ),
        
        shape = guide_legend(
          order = 2,
          title = "Period",
          override.aes = list(
            colour = "black",
            size = 2
          )
        )
      )
  }
  
  
  return(p)
}


# ============================================================
# 10. Four panels
# ============================================================


# ------------------------------------------------------------
# a. Antarctic | Color type
# ------------------------------------------------------------

p_ant_color <- make_trajectory_panel(
  
  pole_use = "Antarctic",
  
  group_type_use = "Color type",
  
  panel_label = "a  Antarctic | Color type",
  
  x_limits = ant_x_limits,
  
  y_limits = c(
    -25,
    220
  ),
  
  show_legend = FALSE,
  
  show_x_title = FALSE,
  
  show_y_title = TRUE
)


# ------------------------------------------------------------
# b. Antarctic | Growth form
# ------------------------------------------------------------

p_ant_growth <- make_trajectory_panel(
  
  pole_use = "Antarctic",
  
  group_type_use = "Growth form",
  
  panel_label = "b  Antarctic | Growth form",
  
  x_limits = ant_x_limits,
  
  y_limits = c(
    -25,
    220
  ),
  
  # Only this panel carries the complete legend.
  show_legend = TRUE,
  
  show_x_title = FALSE,
  
  show_y_title = FALSE
)


# ------------------------------------------------------------
# c. Arctic | Color type
# ------------------------------------------------------------

p_arc_color <- make_trajectory_panel(
  
  pole_use = "Arctic",
  
  group_type_use = "Color type",
  
  panel_label = "c  Arctic | Color type",
  
  x_limits = arc_x_limits,
  
  y_limits = c(
    -50,
    25
  ),
  
  show_legend = FALSE,
  
  show_x_title = TRUE,
  
  show_y_title = TRUE
)


# ------------------------------------------------------------
# d. Arctic | Growth form
# ------------------------------------------------------------

p_arc_growth <- make_trajectory_panel(
  
  pole_use = "Arctic",
  
  group_type_use = "Growth form",
  
  panel_label = "d  Arctic | Growth form",
  
  x_limits = arc_x_limits,
  
  y_limits = c(
    -50,
    25
  ),
  
  show_legend = FALSE,
  
  show_x_title = TRUE,
  
  show_y_title = FALSE
)


# ============================================================
# 11. Combine 2 x 2
#
# Compatible with older patchwork:
# no "& theme(...)" syntax
# ============================================================

top_row <- p_ant_color | p_ant_growth

bottom_row <- p_arc_color | p_arc_growth


p_traj_4panel <- (
  top_row /
    bottom_row
) +
  
  patchwork::plot_layout(
    heights = c(
      1,
      1
    )
  ) +
  
  patchwork::plot_annotation(
    
    title = paste0(
      "Multidimensional redistribution of ",
      "climatically suitable habitats"
    ),
    
    subtitle = paste0(
      "Colours indicate SSP scenarios; circles and triangles represent ",
      "the 2050s and 2090s, respectively. Lines connect each functional ",
      "group from the 2050s to the 2090s under the same scenario. ",
      "Labels identify the 2090s SSP5-8.5 endpoints."
    ),
    
    theme = theme(
      
      plot.title = element_text(
        family = base_family,
        size = base_font_pt + 2,
        face = "bold",
        hjust = 0
      ),
      
      plot.subtitle = element_text(
        family = base_family,
        size = base_font_pt,
        hjust = 0
      )
    )
  )


# ============================================================
# 12. Display
# ============================================================

print(
  p_traj_4panel
)


# ============================================================
# 13. Save
# ============================================================

out_png <- file.path(
  trajectory_figure_dir,
  "Fig6c_multidimensional_redistribution_4panel.png"
)

out_pdf <- file.path(
  trajectory_figure_dir,
  "Fig6c_multidimensional_redistribution_4panel.pdf"
)

out_tiff <- file.path(
  trajectory_figure_dir,
  "Fig6c_multidimensional_redistribution_4panel.tiff"
)


ggsave(
  filename = out_png,
  plot = p_traj_4panel,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = 600
)


ggsave(
  filename = out_pdf,
  plot = p_traj_4panel,
  device = grDevices::cairo_pdf,
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm"
)


ggsave(
  filename = out_tiff,
  plot = p_traj_4panel,
  device = "tiff",
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)


cat(
  "\n========================================\n",
  "Trajectory figure completed.\n",
  "PNG:  ", out_png, "\n",
  "PDF:  ", out_pdf, "\n",
  "TIFF: ", out_tiff, "\n",
  "========================================\n",
  sep = ""
)

