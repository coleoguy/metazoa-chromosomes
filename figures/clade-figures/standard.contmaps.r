## Load libraries
library(ape)
library(diversitree)
library(phangorn)
library(phytools)

###### Parameters ######
# Adjust this value to move the legend up or down.
# More negative values move the legend down, less negative (or positive) moves it up.
# Default is -0.05. Try -0.10 to move down, -0.02 to move up.
legend_y_multiplier <- -0.1

# Directory paths
tree_dir   <- "../../data/trees/"
chrom_dir  <- "../../data/chrome/"
output_dir <- "contmaps/"

# List of plant clades (used to select tree adjustment method)
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", 
            "liliaceae", "solanaceae", "rubiaceae", "passifloraceae")

###### Functions ######

# Function to read and return a tree object for a given clade name
GetTree <- function(x) {
  new_path <- paste0(tree_dir, x, ".new")
  nex_path <- paste0(tree_dir, x, ".nex")
  
  if (file.exists(new_path)) {
    tree <- read.tree(new_path)
  } else if (file.exists(nex_path)) {
    tree <- read.nexus(nex_path)
  } else {
    return(NULL)
  }
  
  tree <- GetCleanTree(tree)
  return(tree)
}

# Function to remove duplicated tip labels from a phylo or multiPhylo object
GetCleanTree <- function(tree) {
  if (inherits(tree, "phylo")) {
    duptips <- which(duplicated(tree$tip.label))
    if (length(duptips) > 0) {
      tree <- drop.tip(tree, duptips)
    }
  } else if (inherits(tree, "multiPhylo")) {
    for (i in seq_along(tree)) {
      duptips <- which(duplicated(tree[[i]]$tip.label))
      if (length(duptips) > 0) {
        tree[[i]] <- drop.tip(tree[[i]], duptips)
      }
    }
  }
  return(tree)
}

# Clean haploid values (handles ranges and comma-separated values)
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  
  for (i in seq_len(nrow(df))) {
    x <- df$haploid[i]
    if (is.numeric(x)) {
      if (x != round(x)) {
        df$haploid[i] <- stoch_round(x)
      }
    } else {
      # Handle dash-separated ranges
      nums <- as.numeric(strsplit(x, split = "-", fixed = TRUE)[[1]])
      if (length(nums) > 2) {
        stop(paste("Looks like you have too many dashes! Check row:", i))
      }
      if (length(nums) == 2) {
        df$haploid[i] <- sample(nums[1]:nums[2], 1)
      } else {
        # Handle comma-separated values
        nums <- as.numeric(strsplit(x, split = ",", fixed = TRUE)[[1]])
        if (length(nums) > 1) {
          df$haploid[i] <- sample(nums, 1)
        }
      }
    }
  }
  df$haploid <- stoch_round(as.numeric(df$haploid))
  return(df)
}

# Makes the tree ultrametric and normalizes branch lengths
AdjustTree <- function(tree, clade) {
  if (!is.ultrametric(tree)) {
    if (clade %in% plants) {
      tree <- nnls.tree(cophenetic(tree), tree, rooted = TRUE)
    } else {
      tree <- chronos(tree)
    }
  }
  tree$edge.length <- tree$edge.length / max(branching.times(tree))
  tree$edge.length[tree$edge.length < 0] <- 1e-9
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
  }
  return(tree)
}

# Determine line width based on number of tips
GetLineWidth <- function(N) {
  if (N < 50) {
    return(4)
  } else if (N < 100) {
    return(3)
  } else if (N < 300) {
    return(1)
  } else {
    return(0.5)
  }
}

###### End Functions ######

# Get all clades from tree files
tree_files <- list.files(tree_dir, full.names = FALSE)
clades <- tools::file_path_sans_ext(basename(tree_files))

for (i in seq_along(clades)) {
  clade <- clades[i]
  output_file <- paste0(output_dir, clade, "_CM.png")
  
  # Skip if output file already exists
  if (file.exists(output_file)) {
    print(paste("Skipping", clade, "- file already exists"))
    next
  }
  
  print(paste("Working on", clade))
  
  # Read chromosome data
  dat <- read.csv(paste0(chrom_dir, clade, ".csv"))
  dat <- CleanHaploid(dat)
  
  # Read tree
  tree <- GetTree(clade)
  if (is.null(tree)) {
    print(paste("No tree found for", clade, "- skipping"))
    next
  }
  
  # Handle multiPhylo objects
  
  if (inherits(tree, "multiPhylo")) {
    tree <- tree[[1]]
  }
  class(tree) <- "phylo"
  
  # Prune tree and data to matching species
  matched <- intersect(tree$tip.label, dat$species)
  tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
  dat <- dat[dat$species %in% matched, c("species", "haploid")]
  dat <- dat[!duplicated(dat$species), ]
  
  # Make tree ultrametric
  tree <- AdjustTree(tree, clade)
  
  # Prepare trait values
  vals <- dat$haploid
  names(vals) <- dat$species
  
  # Create contMap object
  obj <- contMap(tree, vals, plot = FALSE)
  obj <- setMap(obj, hcl.colors(n = 100, palette = "viridis"))
  
  # Plot setup
  N <- length(obj$tree$tip.label)
  lwd <- GetLineWidth(N)
  
  # Save plot
  png(output_file, width = 5, height = 5, units = "in", res = 600)
  plot(obj$tree,
       colors = obj$cols,
       ylim = c(-0.15 * N, N),
       ftype = "off",
       outline = FALSE,
       lwd = lwd)
  
  add.color.bar(leg = 0.5 * max(nodeHeights(obj$tree)),
                cols = obj$cols,
                title = "Haploid Count",
                lims = obj$lims,
                prompt = FALSE,
                x = 0,
                y = legend_y_multiplier * N,
                subtitle = "")
  dev.off()
}