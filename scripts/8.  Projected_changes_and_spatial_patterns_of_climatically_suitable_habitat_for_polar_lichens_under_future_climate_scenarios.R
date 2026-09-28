############################################################
# Figure 5: 2090s SSP5-8.5 only
# Glacier-masked polar lichen suitability rasters
#
# Outputs
# 1) Summary figure with three panels:
#    a) Gain-Loss dumbbell plot based on weighted suitable area (WSA)
#    b) Net-change lollipop plot
#    c) Threshold-defined habitat-loss cell counts
#
# 2) Spatial figure for each functional group:
#    left  = Net suitability change (P_future - P_current)
#    right = Habitat-loss distribution
#
# Definitions
# WSA = sum(P_i * A_i)
# Gain = sum(max(P_future - P_current, 0) * A_i)
# Loss = sum(max(P_current - P_future, 0) * A_i)
# Net  = Gain - Loss
# Percentages are relative to current WSA.
#
# Habitat-loss cell:
# P_current > 0.6 AND P_future < 0.4
#
# Spatial maps:
# - Antarctic study-area outline only
# - Arctic study-area outline + local world-country boundaries
# - Country boundaries are lighter than the study-area outline
# - Habitat-loss cells are drawn as the top layer
############################################################

# ============================================================
# 0. Packages
# ============================================================

install.packages("rnaturalearth")
library(terra)
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)
library(rnaturalearth)



# ============================================================
# 1. Paths and global settings
# ============================================================

masked_root <- "F:/1 MY PROJECTS/Climate Change/RESULT_glacier_masked"
work_root   <- "F:/1 MY PROJECTS/Climate Change/分析过程"

study_area_files <- c(
  "Antarctic" = "F:/1 MY PROJECTS/Climate Change/Study area/antarctica.shp",
  "Arctic"    = "F:/1 MY PROJECTS/Climate Change/Study area/Arctic.shp"
)

# Output folders
fig5_root <- file.path(work_root, "_Fig5_2090s_SSP585")
table_dir <- file.path(fig5_root, "tables")
summary_figure_dir <- file.path(fig5_root, "summary_figure")
spatial_raster_dir <- file.path(fig5_root, "spatial_rasters")
spatial_figure_dir <- file.path(fig5_root, "spatial_figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_figure_dir, recursive = TRUE, showWarnings = FALSE)

# Input suitability settings
probability_divisor <- 1000
alignment_method <- "bilinear"

# Habitat-loss thresholds
current_suitable_threshold <- 0.6
future_unsuitable_threshold <- 0.4

# Figure settings
base_font_pt <- 7
base_family <- "Arial"
map_font_family <- "Arial"

summary_width_mm <- 180
summary_height_mm <- 145
map_width_mm <- 180
map_height_mm <- 95

# Spatial display settings
map_max_plot_cells <- 160000
loss_map_max_plot_cells <- 500000
map_colour_quantile <- 0.995

map_crs <- c(
  "Antarctic" = "EPSG:3031",
  "Arctic"    = "EPSG:3413"
)

terraOptions(progress = 1, memfrac = 0.75)

# ============================================================
# 2. Functional-group order and target scenario
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

# Only one future combination is used in this script.
target_period <- "2090s"
target_scenario <- "SSP585"
target_scenario_label <- "SSP5-8.5"
target_future_name <- "2090_SSP585.img"

# ============================================================
# 3. Helper functions: raster calculations
# ============================================================

raster_range <- function(r) {
  x <- global(r, c("min", "max"), na.rm = TRUE)
  c(
    min = as.numeric(x[1, "min"]),
    max = as.numeric(x[1, "max"])
  )
}

rescale_to_01 <- function(r, file_name) {
  raw_rng <- raster_range(r)
  
  if (!all(is.finite(raw_rng))) {
    stop("Raster contains no finite values: ", file_name)
  }
  
  p <- if (raw_rng["max"] > 1.5) {
    r / probability_divisor
  } else {
    r
  }
  
  scaled_rng <- raster_range(p)
  
  if (
    scaled_rng["min"] < -1e-6 ||
    scaled_rng["max"] > 1 + 1e-6
  ) {
    stop(
      "Suitability values are outside 0-1 after rescaling.\n",
      "File: ", file_name, "\n",
      "Raw range: ", raw_rng["min"], " to ", raw_rng["max"], "\n",
      "Scaled range: ", scaled_rng["min"], " to ", scaled_rng["max"]
    )
  }
  
  clamp(p, lower = 0, upper = 1, values = TRUE)
}

align_to_template <- function(template_raster, input_raster, method = alignment_method) {
  if (isTRUE(compareGeom(template_raster, input_raster, stopOnError = FALSE))) {
    return(input_raster)
  }
  
  if (isTRUE(same.crs(template_raster, input_raster))) {
    resample(input_raster, template_raster, method = method)
  } else {
    project(input_raster, template_raster, method = method, mask = FALSE)
  }
}

g_sum <- function(r) {
  out <- global(r, "sum", na.rm = TRUE)[1, 1]
  if (length(out) == 0 || is.na(out) || !is.finite(out)) return(0)
  as.numeric(out)
}

write_float_raster <- function(r, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  writeRaster(
    r,
    filename,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      NAflag = -9999,
      gdal = c("COMPRESS=LZW", "TILED=YES")
    )
  )
  
  normalizePath(filename, winslash = "/", mustWork = TRUE)
}

write_class_raster <- function(r, filename) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  
  writeRaster(
    r,
    filename,
    overwrite = TRUE,
    wopt = list(
      datatype = "INT1U",
      NAflag = 255,
      gdal = c("COMPRESS=LZW", "TILED=YES")
    )
  )
  
  normalizePath(filename, winslash = "/", mustWork = TRUE)
}

safe_upper_limit <- function(x, probability = map_colour_quantile, fallback = 0.01) {
  x <- x[is.finite(x) & x > 0]
  
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

sample_raster_values <- function(
    raster_file,
    absolute = FALSE,
    sample_size = 100000
) {
  r <- rast(raster_file)
  n_sample <- min(sample_size, ncell(r))
  
  vals <- spatSample(
    r,
    size = n_sample,
    method = "regular",
    na.rm = TRUE,
    values = TRUE,
    as.df = FALSE
  )
  
  vals <- as.numeric(vals)
  vals <- vals[is.finite(vals)]
  
  if (absolute) vals <- abs(vals)
  
  vals
}

# ============================================================
# 4. Helper functions: vector data and plotting projection
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
  
  coords <- as.data.frame(st_coordinates(sf_obj))
  
  level_cols <- grep("^L[0-9]+$", names(coords), value = TRUE)
  
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

load_study_area <- function(shp_file, target_crs) {
  x <- st_read(shp_file, quiet = TRUE)
  x <- st_make_valid(x)
  x <- st_transform(x, target_crs)
  x
}

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
    fact <- ceiling(sqrt(ncell(r2) / max_cells))
    
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

project_binary_df_for_plot <- function(
    r,
    pole,
    max_cells = loss_map_max_plot_cells
) {
  # Explicit 0/1 binary raster before projection.
  r_binary <- ifel(!is.na(r) & r > 0, 1, 0)
  
  r_projected <- project(
    r_binary,
    map_crs[[pole]],
    method = "near",
    mask = FALSE
  )
  
  # max aggregation retains sparse positive cells.
  if (ncell(r_projected) > max_cells) {
    fact <- ceiling(sqrt(ncell(r_projected) / max_cells))
    
    r_projected <- aggregate(
      r_projected,
      fact = fact,
      fun = "max",
      na.rm = TRUE
    )
  }
  
  r_projected <- ifel(r_projected > 0.5, 1, NA)
  names(r_projected) <- "value"
  
  as.data.frame(
    r_projected,
    xy = TRUE,
    na.rm = TRUE
  )
}

# ============================================================
# 5. Check input files
# ============================================================

vector_files_to_check <- study_area_files

if (any(!file.exists(vector_files_to_check))) {
  missing_vector <- vector_files_to_check[!file.exists(vector_files_to_check)]
  stop(
    "Missing vector file(s):\n",
    paste0("  - ", missing_vector, collapse = "\n")
  )
}

required_rasters <- unique(
  unlist(
    lapply(seq_len(nrow(groups)), function(i) {
      c(
        file.path(groups$group_dir[i], "current.img"),
        file.path(groups$group_dir[i], target_future_name)
      )
    })
  )
)

missing_rasters <- required_rasters[!file.exists(required_rasters)]

if (length(missing_rasters) > 0) {
  stop(
    "Missing required suitability raster(s):\n",
    paste0("  - ", missing_rasters, collapse = "\n")
  )
}

# ============================================================
# 6. Load study-area outlines and Arctic geographic boundaries
# ============================================================

study_area_sf <- lapply(
  names(study_area_files),
  function(pole_name) {
    
    load_study_area(
      study_area_files[[pole_name]],
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
# Arctic country boundaries and coastlines from Natural Earth
# ------------------------------------------------------------

cat(
  "\nLoading Natural Earth Arctic boundaries...\n"
)


# 1. International land boundaries
arctic_admin0_sf <- rnaturalearth::ne_download(
  scale = 10,
  type = "admin_0_boundary_lines_land",
  category = "cultural",
  returnclass = "sf"
)


# 2. Coastlines
arctic_coastline_sf <- rnaturalearth::ne_download(
  scale = 10,
  type = "coastline",
  category = "physical",
  returnclass = "sf"
)


# ------------------------------------------------------------
# Transform both directly to Arctic polar stereographic CRS
# ------------------------------------------------------------

arctic_admin0_sf <- sf::st_transform(
  arctic_admin0_sf,
  map_crs[["Arctic"]]
)

arctic_coastline_sf <- sf::st_transform(
  arctic_coastline_sf,
  map_crs[["Arctic"]]
)


# ------------------------------------------------------------
# Crop to slightly larger than Arctic study area
# ------------------------------------------------------------

arctic_bbox <- sf::st_bbox(
  study_area_sf[["Arctic"]]
)


# Add 2% margin to avoid clipping lines exactly on map edge
x_pad <- 0.02 *
  (
    arctic_bbox["xmax"] -
      arctic_bbox["xmin"]
  )

y_pad <- 0.02 *
  (
    arctic_bbox["ymax"] -
      arctic_bbox["ymin"]
  )


arctic_bbox_expanded <- arctic_bbox

arctic_bbox_expanded["xmin"] <-
  arctic_bbox["xmin"] - x_pad

arctic_bbox_expanded["xmax"] <-
  arctic_bbox["xmax"] + x_pad

arctic_bbox_expanded["ymin"] <-
  arctic_bbox["ymin"] - y_pad

arctic_bbox_expanded["ymax"] <-
  arctic_bbox["ymax"] + y_pad


arctic_admin0_sf <- suppressWarnings(
  sf::st_crop(
    arctic_admin0_sf,
    arctic_bbox_expanded
  )
)

arctic_coastline_sf <- suppressWarnings(
  sf::st_crop(
    arctic_coastline_sf,
    arctic_bbox_expanded
  )
)


arctic_admin0_sf <- arctic_admin0_sf[
  !sf::st_is_empty(arctic_admin0_sf),
]

arctic_coastline_sf <- arctic_coastline_sf[
  !sf::st_is_empty(arctic_coastline_sf),
]


cat(
  "  Admin-0 boundary features: ",
  nrow(arctic_admin0_sf),
  "\n",
  sep = ""
)

cat(
  "  Coastline features: ",
  nrow(arctic_coastline_sf),
  "\n",
  sep = ""
)


# ============================================================
# 7. Main calculations: 2090s SSP5-8.5 only
# ============================================================

result_list <- vector("list", nrow(groups))

for (i in seq_len(nrow(groups))) {
  
  pole_i <- groups$pole[i]
  group_type_i <- groups$group_type[i]
  group_name_i <- groups$group_name[i]
  group_dir_i <- groups$group_dir[i]
  out_index_i <- groups$output_index[i]
  
  cat("\n==================================================\n")
  cat(
    "Processing: ",
    pole_i, " | ",
    group_type_i, " | ",
    group_name_i,
    " | 2090s SSP5-8.5\n",
    sep = ""
  )
  
  current_file <- file.path(group_dir_i, "current.img")
  future_file <- file.path(group_dir_i, target_future_name)
  
  current_raw <- rast(current_file)
  future_raw <- rast(future_file)
  
  # Future raster is aligned to the current raster geometry.
  future_raw <- align_to_template(
    current_raw,
    future_raw,
    method = alignment_method
  )
  
  current_p <- rescale_to_01(current_raw, current_file)
  future_p <- rescale_to_01(future_raw, future_file)
  
  # Compare only cells valid in both current and future rasters.
  common_valid <- ifel(
    !is.na(current_p) & !is.na(future_p),
    1,
    NA
  )
  
  current_common <- mask(current_p, common_valid)
  future_common <- mask(future_p, common_valid)
  
  # Actual cell area in km2.
  area_km2 <- cellSize(current_p, unit = "km")
  area_common <- mask(area_km2, common_valid)
  
  # ----------------------------------------------------------
  # WSA gain-loss-net
  # ----------------------------------------------------------
  
  WSA_current <- g_sum(current_common * area_common)
  WSA_future <- g_sum(future_common * area_common)
  
  delta_p <- future_common - current_common
  
  gain_p <- ifel(delta_p > 0, delta_p, 0)
  loss_p <- ifel(delta_p < 0, -delta_p, 0)
  
  Gain <- g_sum(gain_p * area_common)
  Loss <- g_sum(loss_p * area_common)
  Net <- Gain - Loss
  
  gain_percent <- if (WSA_current > 0) {
    100 * Gain / WSA_current
  } else {
    NA_real_
  }
  
  loss_percent <- if (WSA_current > 0) {
    100 * Loss / WSA_current
  } else {
    NA_real_
  }
  
  net_percent <- if (WSA_current > 0) {
    100 * Net / WSA_current
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------
  # Threshold-defined habitat loss
  # ----------------------------------------------------------
  
  current_suitable <- ifel(
    current_common > current_suitable_threshold,
    1,
    NA
  )
  
  habitat_loss <- ifel(
    current_common > current_suitable_threshold &
      future_common < future_unsuitable_threshold,
    1,
    NA
  )
  
  current_suitable_cells <- g_sum(current_suitable)
  habitat_loss_cells <- g_sum(habitat_loss)
  
  current_suitable_area_km2 <- g_sum(current_suitable * area_common)
  habitat_loss_area_km2 <- g_sum(habitat_loss * area_common)
  
  habitat_loss_percent_of_current_area <- if (current_suitable_area_km2 > 0) {
    100 * habitat_loss_area_km2 / current_suitable_area_km2
  } else {
    NA_real_
  }
  
  # Map display raster:
  # 1 = current suitable habitat
  # 2 = threshold-defined habitat-loss cells
  loss_display <- ifel(!is.na(current_suitable), 1, NA)
  loss_display <- ifel(!is.na(habitat_loss), 2, loss_display)
  
  # ----------------------------------------------------------
  # Save intermediate spatial rasters
  # ----------------------------------------------------------
  
  group_out_dir <- file.path(
    spatial_raster_dir,
    pole_i,
    group_type_i,
    group_name_i
  )
  dir.create(group_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  net_file <- write_float_raster(
    delta_p,
    file.path(group_out_dir, "2090_SSP585_Net_deltaP.tif")
  )
  
  loss_display_file <- write_class_raster(
    loss_display,
    file.path(group_out_dir, "2090_SSP585_HabitatLoss_display.tif")
  )
  
  result_list[[i]] <- tibble(
    output_index = out_index_i,
    pole = pole_i,
    group_type = group_type_i,
    group_name = group_name_i,
    period = target_period,
    scenario = target_scenario,
    scenario_label = target_scenario_label,
    current_file = current_file,
    future_file = future_file,
    WSA_current = WSA_current,
    WSA_future = WSA_future,
    Gain = Gain,
    Loss = Loss,
    Net = Net,
    gain_percent = gain_percent,
    loss_percent = loss_percent,
    net_percent = net_percent,
    current_suitable_cells = current_suitable_cells,
    habitat_loss_cells = habitat_loss_cells,
    current_suitable_area_km2 = current_suitable_area_km2,
    habitat_loss_area_km2 = habitat_loss_area_km2,
    habitat_loss_percent_of_current_area = habitat_loss_percent_of_current_area,
    net_file = net_file,
    loss_display_file = loss_display_file
  )
  
  cat(
    "  Gain = ", round(gain_percent, 2), "% | ",
    "Loss = ", round(loss_percent, 2), "% | ",
    "Net = ", round(net_percent, 2), "% | ",
    "Habitat-loss cells = ", habitat_loss_cells,
    "\n",
    sep = ""
  )
}

summary_df <- bind_rows(result_list)

summary_csv <- file.path(
  table_dir,
  "Fig5_2090s_SSP585_summary.csv"
)

write_csv(summary_df, summary_csv, na = "")

# ============================================================
# 8. Summary figure
#    1) Gain-Loss dumbbell
#    2) Net lollipop
#    3) Threshold-defined habitat-loss cell counts
# ============================================================

plot_order_df <- tibble(
  pole = c(rep("Antarctic", 6), rep("Arctic", 6)),
  group_type = c(
    "Color type", "Color type", "Color type",
    "Growth form", "Growth form", "Growth form",
    "Color type", "Color type", "Color type",
    "Growth form", "Growth form", "Growth form"
  ),
  group_name = c(
    "Pale", "Bright", "Dark",
    "Crustose", "Foliose", "Fruticose",
    "Pale", "Bright", "Dark",
    "Crustose", "Foliose", "Fruticose"
  ),
  y = 12:1
)

fig5_df <- plot_order_df %>%
  left_join(
    summary_df %>%
      select(
        pole,
        group_type,
        group_name,
        gain_percent,
        loss_percent,
        net_percent,
        habitat_loss_cells
      ),
    by = c("pole", "group_type", "group_name")
  ) %>%
  mutate(
    gain_percent = as.numeric(gain_percent),
    loss_percent = as.numeric(loss_percent),
    net_percent = as.numeric(net_percent),
    habitat_loss_cells = as.numeric(habitat_loss_cells)
  )

if (
  any(is.na(fig5_df$gain_percent)) ||
  any(is.na(fig5_df$loss_percent)) ||
  any(is.na(fig5_df$net_percent)) ||
  any(is.na(fig5_df$habitat_loss_cells))
) {
  stop("Missing values were found while building the Fig. 5 summary figure.")
}

left_xmin <- -ceiling(max(fig5_df$loss_percent, na.rm = TRUE) / 5) * 5 - 5
left_xmax <- ceiling(max(fig5_df$gain_percent, na.rm = TRUE) / 10) * 10 + 5
net_xmax <- ceiling(max(fig5_df$net_percent, na.rm = TRUE) / 10) * 10 + 5

left_break_min <- floor(left_xmin / 20) * 20
left_break_max <- ceiling(left_xmax / 20) * 20
net_break_max <- ceiling(net_xmax / 20) * 20

# Common theme with visible axes.
theme_fig5 <- theme_classic(
  base_size = base_font_pt,
  base_family = base_family
) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(
      size = base_font_pt,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = base_font_pt,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = base_font_pt,
      face = "bold",
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.35
    ),
    axis.ticks.length = unit(1.2, "mm"),
    plot.title = element_text(
      size = base_font_pt + 0.5,
      face = "bold",
      hjust = 0.5
    ),
    plot.margin = margin(4, 5, 4, 5)
  )

# ------------------------------------------------------------
# Panel 1: Gain-Loss dumbbell
# ------------------------------------------------------------

p_gain_loss <- ggplot(fig5_df, aes(y = y)) +
  geom_hline(
    yintercept = 6.5,
    colour = "grey82",
    linewidth = 0.3
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey55",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      x = -loss_percent,
      xend = gain_percent,
      yend = y
    ),
    colour = "grey45",
    linewidth = 0.45,
    lineend = "round"
  ) +
  geom_point(
    aes(x = -loss_percent),
    shape = 16,
    size = 2.0,
    colour = "#2166AC"
  ) +
  geom_point(
    aes(x = gain_percent),
    shape = 16,
    size = 2.0,
    colour = "#E66101"
  ) +
  geom_text(
    aes(
      x = -loss_percent,
      label = paste0("-", round(loss_percent))
    ),
    hjust = 1.15,
    size = base_font_pt / 2.8,
    colour = "black",
    fontface = "bold"
  ) +
  geom_text(
    aes(
      x = gain_percent,
      label = paste0("+", round(gain_percent))
    ),
    hjust = -0.10,
    size = base_font_pt / 2.8,
    colour = "#E66101",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = left_xmin + 2,
    y = 3.8,
    label = "Loss (%)",
    colour = "#2166AC",
    fontface = "bold",
    hjust = 0,
    size = base_font_pt / 2.2
  ) +
  annotate(
    "text",
    x = max(10, left_xmax * 0.42),
    y = 3.8,
    label = "Gain (%)",
    colour = "#E66101",
    fontface = "bold",
    hjust = 0,
    size = base_font_pt / 2.2
  ) +
  annotate(
    "text",
    x = left_xmin - 6,
    y = 9.5,
    label = "Antarctic",
    angle = 90,
    colour = "#2C5AA0",
    fontface = "bold",
    size = base_font_pt / 2.2
  ) +
  annotate(
    "text",
    x = left_xmin - 6,
    y = 3.5,
    label = "Arctic",
    angle = 90,
    colour = "#1F9BB6",
    fontface = "bold",
    size = base_font_pt / 2.2
  ) +
  scale_y_continuous(
    breaks = fig5_df$y,
    labels = fig5_df$group_name,
    limits = c(0.4, 12.6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(left_xmin, left_xmax),
    breaks = seq(left_break_min, left_break_max, by = 20),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Gain-Loss comparison",
    x = "Area relative to current WSA (%)"
  ) +
  theme_fig5 +
  theme(
    plot.margin = margin(4, 5, 4, 22)
  )

# ------------------------------------------------------------
# Panel 2: Net lollipop
# ------------------------------------------------------------

p_net <- ggplot(fig5_df, aes(y = y)) +
  geom_hline(
    yintercept = 6.5,
    colour = "grey82",
    linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = net_percent,
      yend = y
    ),
    colour = "grey35",
    linewidth = 0.45,
    lineend = "round"
  ) +
  geom_point(
    aes(x = net_percent),
    shape = 16,
    size = 2.0,
    colour = "black"
  ) +
  geom_text(
    aes(
      x = net_percent,
      label = ifelse(
        net_percent >= 0,
        paste0("+", round(net_percent)),
        as.character(round(net_percent))
      )
    ),
    hjust = -0.10,
    size = base_font_pt / 2.8,
    colour = "black",
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = net_xmax * 0.45,
    y = 3.8,
    label = "Net (%)",
    colour = "black",
    fontface = "bold",
    hjust = 0,
    size = base_font_pt / 2.2
  ) +
  scale_y_continuous(
    breaks = fig5_df$y,
    labels = NULL,
    limits = c(0.4, 12.6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0, net_xmax),
    breaks = seq(0, net_break_max, by = 20),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Net change",
    x = "Net change (%)"
  ) +
  theme_fig5 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# ------------------------------------------------------------
# Panel 3: threshold-defined habitat-loss cell counts
# P_current > 0.6 -> P_future < 0.4
# pseudo-log axis allows Antarctic and Arctic values to coexist.
# ------------------------------------------------------------

count_breaks <- c(0, 1, 10, 100, 1000, 3000)
max_count <- max(fig5_df$habitat_loss_cells, na.rm = TRUE)
count_breaks <- count_breaks[count_breaks <= max(3000, max_count * 1.1)]

p_loss_cells <- ggplot(fig5_df) +
  geom_hline(
    yintercept = 6.5,
    colour = "grey82",
    linewidth = 0.3
  ) +
  geom_rect(
    data = fig5_df %>% filter(habitat_loss_cells > 0),
    aes(
      xmin = 0,
      xmax = habitat_loss_cells,
      ymin = y - 0.27,
      ymax = y + 0.27
    ),
    fill = "#B2182B",
    colour = NA
  ) +
  geom_text(
    data = fig5_df %>% filter(habitat_loss_cells > 0),
    aes(
      x = habitat_loss_cells,
      y = y,
      label = comma(habitat_loss_cells)
    ),
    hjust = -0.08,
    size = base_font_pt / 2.8,
    colour = "black",
    fontface = "bold"
  ) +
  geom_text(
    data = fig5_df %>% filter(habitat_loss_cells == 0),
    aes(
      x = 0.6,
      y = y,
      label = "0"
    ),
    hjust = 0,
    size = base_font_pt / 2.8,
    colour = "black",
    fontface = "bold"
  ) +
  scale_y_continuous(
    breaks = fig5_df$y,
    labels = NULL,
    limits = c(0.4, 12.6),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    breaks = count_breaks,
    labels = comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = paste0(
      "Cells with suitability decline\n",
      "from Pcurrent > 0.6 to Pfuture < 0.4"
    ),
    x = "Cell count"
  ) +
  theme_fig5 +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

# ------------------------------------------------------------
# Combine and save summary figure
# ------------------------------------------------------------

p_fig5_summary <- patchwork::wrap_plots(
  p_gain_loss,
  p_net,
  p_loss_cells,
  widths = c(1.55, 0.85, 1.05)
) +
  patchwork::plot_annotation(
    title = "Gain-loss-net comparison under SSP5-8.5 in the 2090s",
    subtitle = paste0(
      "Gain, loss and net change are expressed relative to current weighted suitable area (WSA). ",
      "The right panel shows cells in which suitability declined from Pcurrent > 0.6 to Pfuture < 0.4."
    )
  )

summary_tiff <- file.path(
  summary_figure_dir,
  "Fig5_summary_2090s_SSP585.tiff"
)
summary_pdf <- file.path(
  summary_figure_dir,
  "Fig5_summary_2090s_SSP585.pdf"
)
summary_png <- file.path(
  summary_figure_dir,
  "Fig5_summary_2090s_SSP585.png"
)

ggsave(
  filename = summary_tiff,
  plot = p_fig5_summary,
  device = "tiff",
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = summary_pdf,
  plot = p_fig5_summary,
  device = grDevices::cairo_pdf,
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm"
)

ggsave(
  filename = summary_png,
  plot = p_fig5_summary,
  device = "png",
  width = summary_width_mm,
  height = summary_height_mm,
  units = "mm",
  dpi = 600
)

# ============================================================
# 9. Spatial-map functions
# ============================================================

make_net_map <- function(
    raster_file,
    pole,
    panel_title,
    upper_limit,
    outline_df
) {
  
  r <- rast(
    raster_file
  )
  
  df <- project_continuous_df_for_plot(
    r,
    pole
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
  
  
  # ----------------------------------------------------------
  # Arctic geographic reference
  # ----------------------------------------------------------
  
  if (pole == "Arctic") {
    
    # Coastlines:
    # very light, providing complete continental/island outlines
    p <- p +
      geom_sf(
        data = arctic_coastline_sf,
        inherit.aes = FALSE,
        colour = "grey88",
        linewidth = 0.04
      )
    
    
    # International land boundaries:
    # slightly darker than coastline
    p <- p +
      geom_sf(
        data = arctic_admin0_sf,
        inherit.aes = FALSE,
        colour = "grey82",
        linewidth = 0.045
      )
  }
  
  
  p +
    
    # --------------------------------------------------------
  # Study-area outline
  # darker than Natural Earth reference lines
  # --------------------------------------------------------
  
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
    
    # IMPORTANT:
    # geom_sf requires coord_sf rather than coord_equal
    coord_sf(
      crs = sf::st_crs(
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
      limits = c(
        -upper_limit,
        upper_limit
      ),
      oob = squish,
      name = expression(
        Delta * "P"
      )
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
        family = map_font_family,
        size = base_font_pt,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 2
        )
      ),
      
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.10
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        family = map_font_family,
        size = base_font_pt
      ),
      
      legend.text = element_text(
        family = map_font_family,
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

make_loss_map <- function(
    raster_file,
    pole,
    panel_title,
    outline_df
) {
  
  r <- rast(
    raster_file
  )
  
  
  # ----------------------------------------------------------
  # Combined display raster:
  # 1 = current suitable habitat
  # 2 = habitat-loss cells
  # ----------------------------------------------------------
  
  current_r <- ifel(
    r >= 0.5,
    1,
    NA
  )
  
  loss_r <- ifel(
    r >= 1.5,
    1,
    NA
  )
  
  
  # Project separately so sparse loss cells are retained.
  current_df <- project_binary_df_for_plot(
    current_r,
    pole
  )
  
  loss_df <- project_binary_df_for_plot(
    loss_r,
    pole
  )
  
  
  current_df$map_class <- factor(
    rep(
      "Current suitable habitat",
      nrow(current_df)
    ),
    levels = c(
      "Current suitable habitat",
      "Habitat-loss cells"
    )
  )
  
  
  loss_df$map_class <- factor(
    rep(
      "Habitat-loss cells",
      nrow(loss_df)
    ),
    levels = c(
      "Current suitable habitat",
      "Habitat-loss cells"
    )
  )
  
  
  # ----------------------------------------------------------
  # Bottom layer
  # ----------------------------------------------------------
  
  p <- ggplot() +
    
    geom_raster(
      data = current_df,
      aes(
        x = x,
        y = y,
        fill = map_class
      )
    )
  
  
  # ----------------------------------------------------------
  # Arctic geographic reference
  # ----------------------------------------------------------
  
  if (pole == "Arctic") {
    
    # Coastline
    p <- p +
      geom_sf(
        data = arctic_coastline_sf,
        inherit.aes = FALSE,
        colour = "grey88",
        linewidth = 0.04
      )
    
    
    # International boundaries
    p <- p +
      geom_sf(
        data = arctic_admin0_sf,
        inherit.aes = FALSE,
        colour = "grey82",
        linewidth = 0.045
      )
  }
  
  
  p +
    
    # --------------------------------------------------------
  # Study-area outline
  # --------------------------------------------------------
  
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
    
    # --------------------------------------------------------
  # TOP layer:
  # threshold-defined habitat-loss cells
  # --------------------------------------------------------
  
  geom_raster(
    data = loss_df,
    aes(
      x = x,
      y = y,
      fill = map_class
    )
  ) +
    
    coord_sf(
      crs = sf::st_crs(
        map_crs[[pole]]
      ),
      datum = NA,
      expand = FALSE
    ) +
    
    scale_fill_manual(
      values = c(
        "Current suitable habitat" = "grey85",
        "Habitat-loss cells" = "#B2182B"
      ),
      
      breaks = c(
        "Current suitable habitat",
        "Habitat-loss cells"
      ),
      
      limits = c(
        "Current suitable habitat",
        "Habitat-loss cells"
      ),
      
      drop = FALSE,
      
      name = NULL
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
        family = map_font_family,
        size = base_font_pt,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 2
        )
      ),
      
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.10
      ),
      
      legend.position = "bottom",
      
      legend.text = element_text(
        family = map_font_family,
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
# 10. Draw and save spatial figures: 2090s SSP5-8.5 only
# ============================================================

graphics.off()

for (i in seq_len(nrow(groups))) {
  
  pole_i <- groups$pole[i]
  group_type_i <- groups$group_type[i]
  group_name_i <- groups$group_name[i]
  out_index_i <- groups$output_index[i]
  
  cat(
    "\nDrawing spatial figure: ",
    pole_i, " | ",
    group_type_i, " | ",
    group_name_i,
    " | 2090s SSP5-8.5\n",
    sep = ""
  )
  
  row_i <- summary_df %>%
    filter(
      pole == pole_i,
      group_type == group_type_i,
      group_name == group_name_i
    )
  
  if (nrow(row_i) != 1L) {
    stop(
      "Expected exactly one calculated result for: ",
      pole_i, " | ", group_type_i, " | ", group_name_i
    )
  }
  
  outline_df <- outline_df_list[[pole_i]]
  
  net_limit <- safe_upper_limit(
    sample_raster_values(
      row_i$net_file[1],
      absolute = TRUE
    )
  )
  
  net_plot <- make_net_map(
    raster_file = row_i$net_file[1],
    pole = pole_i,
    panel_title = "Net suitability change",
    upper_limit = net_limit,
    outline_df = outline_df
  )
  
  loss_plot <- make_loss_map(
    raster_file = row_i$loss_display_file[1],
    pole = pole_i,
    panel_title = "Habitat-loss distribution",
    outline_df = outline_df
  )
  
  combined_spatial <- patchwork::wrap_plots(
    net_plot,
    loss_plot,
    widths = c(1, 1)
  ) +
    patchwork::plot_annotation(
      title = paste0(
        pole_i, " | ",
        group_type_i, " | ",
        group_name_i, " | 2090s SSP5-8.5"
      ),
      subtitle = paste0(
        "Left: net suitability change. ",
        "Right: habitat-loss distribution ",
        "(grey = current suitable habitat; red = habitat-loss cells)."
      )
    )
  
  out_stem <- paste(
    out_index_i,
    pole_i,
    gsub(" ", "_", group_type_i),
    group_name_i,
    "2090s_SSP585_Net_HabitatLoss",
    sep = "_"
  )
  
  map_tiff <- file.path(
    spatial_figure_dir,
    paste0(out_stem, ".tiff")
  )
  map_pdf <- file.path(
    spatial_figure_dir,
    paste0(out_stem, ".pdf")
  )
  map_png <- file.path(
    spatial_figure_dir,
    paste0(out_stem, ".png")
  )
  
  ggsave(
    filename = map_tiff,
    plot = combined_spatial,
    device = "tiff",
    width = map_width_mm,
    height = map_height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw"
  )
  
  ggsave(
    filename = map_pdf,
    plot = combined_spatial,
    device = grDevices::cairo_pdf,
    width = map_width_mm,
    height = map_height_mm,
    units = "mm"
  )
  
  ggsave(
    filename = map_png,
    plot = combined_spatial,
    device = "png",
    width = map_width_mm,
    height = map_height_mm,
    units = "mm",
    dpi = 600
  )
  
  cat(
    "  Habitat-loss cells: ", row_i$habitat_loss_cells[1], "\n",
    "  Saved: ", out_stem, "\n",
    sep = ""
  )
  
  graphics.off()
}

# ============================================================
# 11. Completion message
# ============================================================

cat("\n==================================================\n")
cat("Completed Figure 5 calculations and plotting.\n")
cat("Scenario: 2090s | SSP5-8.5 only\n")
cat("Summary table:\n  ", summary_csv, "\n", sep = "")
cat("Summary figure:\n  ", summary_png, "\n", sep = "")
cat("Spatial figures directory:\n  ", spatial_figure_dir, "\n", sep = "")
cat("Intermediate spatial rasters:\n  ", spatial_raster_dir, "\n", sep = "")
cat("==================================================\n")


