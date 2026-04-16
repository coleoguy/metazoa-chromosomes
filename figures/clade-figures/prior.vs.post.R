library(ape)
library(coda)
library(diversitree)

# --- Setup ---
mcmc_path  <- "../../results/exponential-prior-full-model/mentor_results/"
files      <- list.files(mcmc_path, full.names = TRUE)
ages       <- read.csv("../../data/clade.ages.csv")
output_dir <- "clade_audit_plots" # Folder to save images
dir.create(output_dir, showWarnings = FALSE)

# --- 1. Generate the "Unit" Prior ONCE ---
# (Same logic as before: Sum of 2 Exponentials)
lik.prior <- function(p) { return(0) }
res.prior <- diversitree::mcmc(lik = lik.prior, x.init = runif(2), 
                  prior = make.prior.exponential(2), 
                  nsteps = 10000, w = 1, lower = 0)
unit_prior <- res.prior[[2]] + res.prior[[3]]

# --- 2. The Plotting Loop ---
for (f in files) {
  # Parse Clade Name
  clade_name <- strsplit(basename(f), split=".", fixed=TRUE)[[1]][1]
  
  # Find Age
  age_idx <- match(tolower(clade_name), tolower(ages$Clade))
  if (is.na(age_idx)) {
    warning(paste("Skipping", clade_name, "- No age found."))
    next
  }
  age <- ages$Age[age_idx]
  
  # Load Posterior (Empirical)
  dat <- read.csv(f)
  dat <- dat[dat$i > 200, ] # Burnin
  post_rate <- (dat$asc1 + dat$desc1) / age
  
  # Scale Prior by Age
  prior_rate <- unit_prior / age
  
  # Calculate Densities
  d_post <- density(post_rate, from=0)
  d_prior <- density(prior_rate, from=0)
  
  # Setup Plot Range (Zoom in on the data, but keep some prior context)
  x_max <- max(quantile(post_rate, 0.99), quantile(prior_rate, 0.80))
  y_max <- max(d_post$y, d_prior$y) * 1.1
  
  # --- SAVE TO PNG ---
  png(filename = paste0(output_dir, "/", clade_name, "_audit.png"),
      width = 3, height = 3, units = "in", res = 300)
  
  par(mar = c(4, 4, 2, 1))
  
  # 1. Empty Plot
  plot(1, type="n", xlim=c(0, x_max), ylim=c(0, y_max),
       xlab = "Dysploidy Rate (events / Myr)", 
       ylab = "Density",
       main = paste("Clade:", tools::toTitleCase(clade_name)),
       bty="l")
  
  # 2. Draw Prior (The Null)
  polygon(c(d_prior$x, 0), c(d_prior$y, 0), 
          col = "gray90", border = "gray60", lty=2)
  
  # 3. Draw Posterior (The Signal)
  # You can customize the color here (e.g., using your class_cols lookup)
  polygon(c(d_post$x, 0), c(d_post$y, 0), 
          col = adjustcolor("#1D3557", alpha.f=0.6), border = "#1D3557", lwd=2)
  
  # 4. Legend
  legend("topright", legend=c("Posterior (Data)", "Prior (Null)"),
         fill=c(adjustcolor("#1D3557", alpha.f=0.6), "gray90"),
         border=c("#1D3557", "gray60"), bty="n")
  
  # 5. Add Stats (Optional but powerful)
  # Show the median rate on the plot
  text(x = x_max, y = y_max * 0.5, 
       labels = paste0("Med. Rate: ", round(median(post_rate), 4)), 
       adj = 1, cex = 0.8)
  
  dev.off()
  
  message(paste("Plotted:", clade_name))
}
