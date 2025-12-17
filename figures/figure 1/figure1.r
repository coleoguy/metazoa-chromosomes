# ============================
# Tree + Haploid Beeswarm + Clade Labels (base R, aligned to tip ends)
# ============================

# ---- 1. Libraries ----
library(ape)
library(dplyr)
library(purrr)
library(stringr)
library(readr)

# ---- 2. Load haploid chromosome data ----
file_list <- list.files(path = "../../data/chrome", pattern = "\\.csv$", full.names = TRUE)

all_data <- list()
for (file in file_list) {
  df <- read.csv(file)
  if (!"haploid" %in% names(df)) {
    cat("Missing 'haploid' column in:", basename(file), "\n")
    next
  }
  file_name <- tools::file_path_sans_ext(basename(file))
  temp_df <- data.frame(clade = tolower(file_name), haploid = df$haploid)
  all_data[[length(all_data) + 1]] <- temp_df
}
dat <- do.call(rbind, all_data)

# ---- 3. Load tree ----
tree <- read.tree("cladetree.new")
tree$tip.label <- tolower(tree$tip.label)

# ---- 4. Clean haploid values ----
dat_clean <- dat %>%
  mutate(
    haploid = str_trim(haploid),
    haploid = map_dbl(haploid, ~ {
      vals <- str_split(.x, ",|–|-")[[1]] |> str_trim()
      nums <- suppressWarnings(as.numeric(vals))
      nums <- nums[!is.na(nums) & nums > 0]
      if (length(nums) == 0) return(NA_real_)
      if (length(nums) == 1) return(nums)
      if (length(nums) == 2 && diff(nums) > 0) return(runif(1, nums[1], nums[2]))
      return(sample(nums, 1))
    })
  ) %>%
  filter(is.finite(haploid), haploid > 0)

# ---- 5. Match data to tree tips ----
final_df <- dat_clean %>%
  filter(clade %in% tree$tip.label) %>%
  mutate(
    clade       = factor(clade, levels = tree$tip.label),
    log_haploid = log10(haploid)
  )

# ---- 5b. Add higher classification from CSV ----
class_df <- read.csv("higher_class.csv", stringsAsFactors = FALSE)
class_df <- class_df %>%
  mutate(Clade = tolower(str_trim(Clade)))

final_df <- final_df %>%
  mutate(clade = tolower(str_trim(as.character(clade)))) %>%
  left_join(class_df, by = c("clade" = "Clade")) %>%
  rename(higher_class = Higher.Classification)

# ---- 6. Calculate species counts for labels ----
counts <- final_df %>%
  count(clade, name = "count") %>%
  mutate(
    label_text = paste0(tools::toTitleCase(as.character(clade)), " (n=", count, ")")
  )

# ---- 7. Color palette for higher classes ----
class_cols <- c(
  "Bryophyta"      = "#1D3557", 
  "Pteridophyta"   = "#2E8B57", 
  "Gymnospermae"   = "#0B4218", 
  "Angiospermae"   = "#84D957", 
  "Fungi"          = "#E6AB02", 
  "Insecta"        = "#8E3B9E", 
  "Arachnida"      = "#D00000", 
  "Reptilia"       = "#F48C06", 
  "Mammalia"       = "#6096BA",  
  "Amphibia"       = "#48CAE4", 
  "Actinopterygii" = "#0077B6", 
  "Chondrichthyes" = "#023E8A" 
)

plot_alpha <- 0.6

# ---- 8. Geometry: get tree coordinates & set panel widths ----
# Invisible plot to get base x-range
plot(tree, direction = "rightwards", cex = 0.6, plot = FALSE)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

x_range    <- range(lp$xx)
bees_width <- diff(x_range) * 1.5        # width of beeswarm panel
label_gap  <- diff(x_range) * 0.1        # gap between beeswarm and labels

# How much space do labels need?
label_margin <- max(strwidth(counts$label_text, cex = 0.6)) * 1.2

# Provisional xlim including extra space to the right
xlim_new <- c(
  x_range[1],
  x_range[2] + bees_width + label_gap + label_margin
)

par(xpd = NA, mar = c(5, 4, 4, 2) + 0.1)

# ---- 9. Plot tree (no tip labels) ----
shrink_factor <- 0.4   # try 0.2–0.6
tree$edge.length <- tree$edge.length * shrink_factor


plot(tree,
     direction      = "rightwards",
     cex            = 0.6,
     x.lim          = xlim_new,
     show.tip.label = FALSE)

# Get coordinates from this final tree plot
lp        <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y     <- lp$yy[1:Ntip(tree)]
x_tip_max <- max(lp$xx[1:Ntip(tree)])   # x-position at ends of tips

# Start beeswarm exactly at tip ends
x_start_bees <- x_tip_max
x_end_bees   <- x_start_bees + bees_width
x_start_lab  <- x_end_bees + label_gap

# ---- 10. Prepare beeswarm coordinates ----
global_min <- min(final_df$log_haploid, na.rm = TRUE)
global_max <- max(final_df$log_haploid, na.rm = TRUE)

scale_x <- function(val, global_min, global_max, start_x, width) {
  norm_val <- (val - global_min) / (global_max - global_min)
  start_x + norm_val * width
}

tip_index <- setNames(1:Ntip(tree), tree$tip.label)

final_df <- final_df %>%
  mutate(
    tip_i = tip_index[as.character(clade)],
    y0    = tip_y[tip_i],
    x_bee = scale_x(log_haploid, global_min, global_max, x_start_bees, bees_width),
    x_bee = jitter(x_bee, amount = bees_width * 0.01),
    y_bee = jitter(y0,    amount = 0.2)
  )

counts <- counts %>%
  mutate(
    tip_i = tip_index[as.character(clade)],
    y0    = tip_y[tip_i]
  )

# ---- 11. Legend order based on vertical clade order (top → bottom) ----
tip_order_bottom_top <- tree$tip.label          # ape: 1 = bottom, N = top
tip_order_top_bottom <- rev(tip_order_bottom_top)

clades_with_data <- unique(as.character(final_df$clade))
higher_map       <- setNames(class_df$Higher.Classification, class_df$Clade)

higher_order <- tip_order_top_bottom[tip_order_top_bottom %in% clades_with_data]
higher_order <- higher_map[higher_order]
higher_order <- unique(higher_order[!is.na(higher_order)])

legend_cols  <- class_cols[higher_order]
legend_fills <- sapply(legend_cols, function(x) adjustcolor(x, alpha.f = plot_alpha))

# ---- 12. Draw haploid baselines aligned to tip ends ----
segments(
  x0  = x_start_bees,
  y0  = tip_y,
  x1  = x_end_bees,
  y1  = tip_y,
  col = "grey80",
  lwd = 0.5
)

# ---- 13. Beeswarm points ----
point_cols <- with(final_df, ifelse(
  !is.na(higher_class) & higher_class %in% names(class_cols),
  class_cols[higher_class],
  "grey50"
))

points(
  x   = final_df$x_bee,
  y   = final_df$y_bee,
  pch = 16,
  cex = 0.4,
  col = adjustcolor(point_cols, alpha.f = 0.1)
)

# ---- 14. Clade labels on the far right ----
text(
  x      = x_start_lab,
  y      = counts$y0,
  labels = counts$label_text,
  adj    = c(0, 0.5),
  cex    = 0.9
)

# ---- 15. X-axis for beeswarm (log10 haploid) ----
axis_y_pos <- 0.2

tick_vals <- log10(c(2, 5, 10, 20, 50, 100, 200))
tick_vals <- tick_vals[tick_vals >= global_min & tick_vals <= global_max]
tick_locs <- scale_x(tick_vals, global_min, global_max, x_start_bees, bees_width)

# Draw the axis line across the *entire* beeswarm panel
segments(
  x0 = x_start_bees,
  y0 = axis_y_pos,
  x1 = x_end_bees,
  y1 = axis_y_pos,
  lwd = 1
)

axis(
  side      = 1,
  at        = tick_locs,
  labels    = c(2, 5, 10, 20, 50, 100, 200)[seq_len(length(tick_vals))],
  pos       = axis_y_pos,
  lwd       = 0,
  lwd.ticks = 1,
  cex.axis  = 1.2
)

mtext(
  "Haploid chromosome number",
  side = 1,
  line = 2,
  at   = x_start_bees + (bees_width / 2),
  cex  = 1.6
)

# ---- 16. Legend ----
legend(
  "bottomleft",
  legend    = higher_order,
  col       = legend_cols,
  pt.bg     = legend_fills,
  pch       = 22,
  pt.cex    = 1.5,
  cex       = 0.85,
  bty       = "n",
  y.intersp = 0.6,
  x.intersp = 0.6,
  inset     = c(-0.06, 0)   # <- push left (try -0.05 to -0.15)
)
