library(phytools)

# --- CONFIGURATION ---
tree_dir  <- "../../data/trees"
chrom_dir <- "../../data/chrome"

# Get all files
files <- list.files(tree_dir, pattern = "\\.(tre|nex|new|nexus)$", full.names = TRUE)

# Pause between plots so you can examine them in the window
par(ask = TRUE)

# --- MAIN LOOP ---
# CHANGED: Using numeric iterator as requested for easier debugging
for (i in 1:length(files)) {
  
  # Set current file based on index
  f <- files[i]
  
  # 1. LOAD TREE
  name <- tools::file_path_sans_ext(basename(f))
  tr   <- tryCatch(read.nexus(f), error = function(e) read.tree(f))
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  
  # 2. MATCH DATA
  csv_path <- file.path(chrom_dir, paste0(name, ".csv"))
  dat  <- read.csv(csv_path)
  keep <- intersect(tr$tip.label, dat$species)
  tr_clean <- keep.tip(tr, keep)
  tr_clean$edge.length[tr_clean$edge.length <= 0] <- 1e-6
  vals <- setNames(dat$haploid[match(keep, dat$species)], keep)
  
  # 3. PREP MAP (Default Palette)
  obj <- contMap(tr_clean, vals, type="fan",plot = T)
  
  # 4. RENDER TO PLOTTING WINDOW
  my_lwd <- max(1.5, 600 / length(keep))
  message(paste("Plotting [", i, "]: ", name, sep=""))
  
  # Layout: 85% Tree, 15% Legend
  layout(matrix(c(1, 2), nrow = 2), heights = c(0.85, 0.15))
  
  # Plot Tree
  par(mar = c(0, 0, 0, 0))
  plot(obj, type = "fan", lwd = my_lwd, outline = FALSE, ftype = "off", legend = FALSE)
  
  # Plot Legend Canvas
  par(mar = c(0, 0, 0, 0))
  plot(0, type = "n", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  
  # Draw Legend (Using your fixed coordinates)
  add.color.bar(leg = 0.5, 
                cols = obj$cols, 
                title = "Haploid Number", 
                lims = obj$lims, 
                digits = 1, 
                prompt = FALSE, 
                x = 0,      # Left-aligned as per your fix
                y = 0.8,    # Lifted up as per your fix
                lwd = 10, 
                fsize = 0.8, 
                subtitle = "")       
}

par(ask = FALSE)
message("Job done.")