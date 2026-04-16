library(ape)
library(coda)

# --- 1. Setup ---
tree       <- read.tree("../figure-1/cladetree.new")
higher_df  <- read.csv("../figure-1/higher_class.csv")
mcmc_files <- list.files("../../results/exponential-prior-full-model/mentor_results/", full.names = TRUE)
ages       <- read.csv("../../data/clade.ages.csv")

# VISUALIZATION SETTINGS
plot_alpha <- 0.4             
min_rate_threshold <- 1e-05   
nbuffer <- 0.4  # <--- TWEAK THIS: Larger number moves the plots further right

hc_lookup  <- setNames(higher_df$Higher.Classification, tolower(higher_df$Clade))

#--- FIX 1: Drop magnolicaceae from tree
tree <- drop.tip(tree, "Magnoliaceae")

# --- FIX 2: Optimized Color Palette (3 Greens + Maximal Divergence) ---
class_cols <- c(
  # Plants (Green Gradient)
  "Bryophyta"      = "#1D3557", 
  "Pteridophyta"   = "#2E8B57", 
  "Gymnospermae"   = "#0B4218", 
  "Angiospermae"   = "#84D957", 
  # Fungi
  "Fungi"          = "#E6AB02", 
  # Invertebrates
  "Insecta"        = "#8E3B9E", 
  "Arachnida"      = "#D00000", 
  # Vertebrates
  "Reptilia"       = "#F48C06", 
  "Mammalia"       = "#6096BA",  
  "Amphibia"       = "#48CAE4", 
  "Actinopterygii" = "#0077B6", 
  "Chondrichthyes" = "#023E8A" 
)

# --- 2. Data Processing ---
dens_data  <- list()
global_rng <- c(log10(min_rate_threshold), -Inf) 

for (i in 1:length(mcmc_files)) {
  clade <- strsplit(basename(mcmc_files[i]), split=".", fixed=T)[[1]][1]
  age_idx <- match(tolower(clade), tolower(ages$Clade))

  age <- ages$Age[age_idx]
  dat <- read.csv(mcmc_files[i])
  dat <- dat[dat$i > 100, c(-1, -ncol(dat))]
  
  vals <- rowSums(dat[,1:2])
  rate_vals <- vals / age
  
  valid_rates <- rate_vals[rate_vals > 0]

  log_rates <- log10(valid_rates)
  
  d <- density(log_rates, cut = 0.5) 
  
  global_rng[2] <- max(global_rng[2], max(d$x))
  
  dens_data[[clade]] <- list(x = d$x, y = d$y, mean = mean(log_rates))
}

# --- 3. Geometry ---
plot(tree, direction = "rightwards", cex = 0.6, plot = FALSE)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

x_range    <- range(lp$xx)
dens_width <- diff(x_range) * 1.5      
offset     <- diff(x_range) * nbuffer # <--- USING YOUR BUFFER VARIABLE HERE
xlim_new   <- c(x_range[1], x_range[2] + offset + dens_width)
tip_y      <- lp$yy[1:Ntip(tree)]

scale_x <- function(val, global_min, global_max, start_x, width) {
  norm_val <- (val - global_min) / (global_max - global_min)
  start_x + norm_val * width
}

# --- 4. Plotting ---
par(xpd = NA, mar = c(5, 4, 4, 2) + 0.1)

# Plot Tree
plot(tree, direction = "rightwards", cex = 0.6, x.lim = xlim_new)

# Pre-calculate text widths for guide lines
label_widths <- strwidth(tree$tip.label, cex = 0.6)
text_buffer  <- diff(x_range) * 0.02 

bump_height  <- 4
x_start_dens <- x_range[2] + offset

# --- Grid Lines ---
log_tick_vals <- pretty(global_rng, n = 5)
log_tick_vals <- log_tick_vals[log_tick_vals >= global_rng[1]]
tick_locs     <- scale_x(log_tick_vals, global_rng[1], global_rng[2], x_start_dens, dens_width)

segments(x0 = tick_locs, y0 = 0, 
         x1 = tick_locs, y1 = Ntip(tree) + 1, 
         col = "lightgray", lty = "dotted", lwd = 0.8)

# REVERSE LOOP
for (i in seq(Ntip(tree), 1, by = -1)) {
  clade <- tree$tip.label[i]
  
  dens <- dens_data[[clade]]
  if (is.null(dens)) {
    matches <- names(dens_data)[tolower(names(dens_data)) == tolower(clade)]
    if (length(matches) > 0) dens <- dens_data[[matches[1]]]
  }
  if (is.null(dens)) next
  
  grp <- hc_lookup[tolower(clade)]
  base_col <- if (!is.na(grp) && grp %in% names(class_cols)) class_cols[grp] else "grey50"
  final_col <- adjustcolor(base_col, alpha.f = plot_alpha)
  
  keep_idx <- dens$x >= global_rng[1]
  x_trunc  <- dens$x[keep_idx]
  y_trunc  <- dens$y[keep_idx]
  
  if (length(x_trunc) < 2) next
  
  y0 <- tip_y[i]
  x_scaled <- scale_x(x_trunc, global_rng[1], global_rng[2], x_start_dens, dens_width)
  y_scaled <- y0 + (y_trunc / max(dens$y)) * bump_height
  
  # Draw Guide Line
  x_line_start <- lp$xx[i] + label_widths[i] + text_buffer
  x_poly_start <- x_scaled[1]
  
  if (x_poly_start > x_line_start) {
    segments(x0 = x_line_start, y0 = y0, 
             x1 = x_poly_start, y1 = y0, 
             col = "gray50", lty = 3, lwd = 1)
  }
  
  # Draw Polygon
  poly_x <- c(x_scaled[1], x_scaled, tail(x_scaled, 1))
  poly_y <- c(y0, y_scaled, y0)
  
  polygon(poly_x, poly_y, border = NA, col = final_col)
  lines(x_scaled, y_scaled, col = base_col, lwd = 1)
  
  # Draw Mean Line
  if (dens$mean >= global_rng[1]) {
    x_mean <- scale_x(dens$mean, global_rng[1], global_rng[2], x_start_dens, dens_width)
    mean_height <- (approx(x_trunc, y_trunc, xout=dens$mean)$y / max(dens$y)) * bump_height
    segments(x_mean, y0, x_mean, y0 + mean_height, col = base_col, lwd = 2)
  }
}

# --- 5. Decoration ---
axis_y_pos <- 0.2

# Ticks
log_tick_vals <- pretty(global_rng, n = 5)
log_tick_vals <- log_tick_vals[log_tick_vals >= global_rng[1]]
tick_locs     <- scale_x(log_tick_vals, global_rng[1], global_rng[2], x_start_dens, dens_width)

# Axis Line
segments(x0 = scale_x(global_rng[1], global_rng[1], global_rng[2], x_start_dens, dens_width), 
         y0 = axis_y_pos, 
         x1 = tail(tick_locs, 1), 
         y1 = axis_y_pos, 
         lwd = 1)

# Smart Labels
smart_labels <- sapply(log_tick_vals, function(x) {
  if (x == 0) {
    return(expression(1))
  } else {
    return(parse(text=paste("10^", as.integer(x), sep="")))
  }
})

# Ticks & Labels
axis(1, at = tick_locs, labels = smart_labels, pos = axis_y_pos, 
     lwd = 0, lwd.ticks = 1, cex.axis = 0.8)

# Title
mtext("Dysploidy Rate (events / Myr)", side = 1, line = 1.5, 
      at = x_start_dens + (dens_width / 2), cex = 0.9)

# Legend
legend_fills <- sapply(class_cols, function(x) adjustcolor(x, alpha.f = plot_alpha))
legend("bottomleft",
       inset = c(-.03,0),
       legend = names(class_cols), 
       col = class_cols,       
       pt.bg = legend_fills,   
       pch = 22, 
       pt.cex = 1.5, 
       bty = "n", 
       cex = 0.65)
