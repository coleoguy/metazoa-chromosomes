## =========================
## 0. Packages
## =========================
library(ape)

## =========================
## 1. Read and plot tree
## =========================

tree      <- read.tree("../figure 1/cladetree.new")

# quick plot to get coordinates
plot(tree, direction = "rightwards", cex = 0.6)
lp      <- get("last_plot.phylo", envir = .PlotPhyloEnv)
x_range <- range(lp$xx)

dens_width <- diff(x_range) * 1.1
offset     <- diff(x_range) * 0.2

xlim_new <- c(x_range[1], x_range[2] + offset + dens_width)

# final tree plot with extra space for densities
plot(tree, direction = "rightwards", cex = 0.6, x.lim = xlim_new)

lp    <- get("last_plot.phylo", envir = .PlotPhyloEnv)
x_max <- max(lp$xx)
tip_y <- lp$yy[1:Ntip(tree)]


## =========================
## 1b. Higher classification + colors
## =========================
higher_df <- read.csv("../figure 1/higher_class.csv")

# assume columns: "Clade" and "Higher Classification"
hc_lookup <- higher_df[["Higher.Classification"]]
names(hc_lookup) <- tolower(higher_df[["Clade"]])

class_cols <- c(
  "Chondrichthyes" = "#414487",
  "Actinopterygii" = "#2A788E",
  "Amphibia"       = "#22A884",
  "Mammalia"       = "#7AD151",
  "Reptilia"       = "#F8961E",
  "Arachnida"      = "#D1495B",
  "Insecta"        = "#8E3B9E",
  "Angiosperm"     = "#1F968B",
  "Gymnosperms"    = "#89C2D9",
  "Pteridophytes"  = "#5E4FA2",
  "Bryophyta"      = "#277DA1"
)


## =========================
## 2. Read all *_mcmc.csv and get log(asc1[i > 500])
## =========================
mcmc_dir   <- file.path("results", "Results_mcmc")
mcmc_files <- list.files("../../results/trainee_results", full.names = TRUE)
# we will use these ages to get rates into units MY
ages <- read.csv("../../data/clade.ages.csv")

# global scale for log(asc1/age)
global_min <- Inf
global_max <- -Inf

for (f in mcmc_files) {
  dat  <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  clade <- strsplit(basename(f), "_")[[1]][1]
  age   <- ages$Age[tolower(ages$Clade) == tolower(clade)]
  
  log_vals <- log(dat$asc1[dat$i > 500 & dat$asc1 > 0] / age)
  
  if (length(log_vals) > 0) {
    global_min <- min(global_min, min(log_vals))
    global_max <- max(global_max, max(log_vals))
  }
}

par(xpd = NA)          # allow drawing beyond plot region
bump_height <- 1.1     # vertical size of each density bump


## =========================
## 3. Draw colored log-densities + mean lines + fill
## =========================

for (f in mcmc_files) {
  dat <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  age <- ages$Age[tolower(ages$Clade) == strsplit(basename(f), "_")[[1]][1]]
  # log-transform in one line
  log_vals <- log(dat$asc1[dat$i > 500 & dat$asc1 > 0]/age)
  if (length(log_vals) < 2) next  # need ≥2 points for density()
  
  clade <- sub("_mcmc\\.csv$", "", basename(f))
  idx   <- match(tolower(clade), tolower(tree$tip.label))
  if (is.na(idx)) next
  
  y0 <- tip_y[idx]
  
  hc_name  <- hc_lookup[tolower(clade)]
  col_line <- class_cols[hc_name]
  if (is.na(col_line)) col_line <- "black"
  
  # density of log(asc1[i > 500])
  d <- density(log_vals)
  
  # scale density x into reserved band on the right
  x_scaled <- x_max + offset +
    (d$x - global_min) / (global_max - global_min) * dens_width
  
  # scale density y to a bump around y0
  y_scaled <- y0 + d$y / max(d$y) * bump_height
  
  # light fill under curve
  col_fill <- adjustcolor(col_line, alpha.f = 0.3)
  polygon(
    x = c(x_scaled[1], x_scaled, x_scaled[length(x_scaled)]),
    y = c(y0,          y_scaled, y0),
    border = NA,
    col    = col_fill
  )
  
  # colored density curve
  lines(x_scaled, y_scaled, col = col_line, lwd = 1.2)
  
  # grey baseline
  segments(x_max + offset, y0,
           x_max + offset + dens_width, y0,
           col = "grey80", lwd = 0.5)
  
  # vertical mean line (in log space)
  mean_log      <- mean(log_vals)
  x_mean_scaled <- x_max + offset +
    (mean_log - global_min) / (global_max - global_min) * dens_width
  
  segments(x_mean_scaled, y0,
           x_mean_scaled, y0 + bump_height,
           col = col_line, lwd = 1)
}


## =========================
## 4. Legend
## =========================
legend("bottomleft",
       legend    = names(class_cols),
       col       = class_cols,
       pch       = 16,
       pt.cex    = 0.7,
       cex       = 0.7,
       y.intersp = 0.6,
       x.intersp = 0.6,
       bty       = "n",
       title     = "Higher classification")

