## Load libraries
library(ape)           
library(diversitree)   
library(phangorn)
library(phytools)

###### Functions ######
# function to read and return a tree object for a given clade name
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0("../../data/trees/", x, ".new")  # <-- change path as needed
  nex_path <- paste0("../../data/trees/", x, ".nex")  # <-- change path as needed

  
  # check which file exists and read accordingly
  if (file.exists(new_path)) {
    tree <- read.tree(new_path)
  } else if (file.exists(nex_path)) {
    tree <- read.nexus(nex_path)
  } else {
    return(NULL)  # skip if no matching tree file found
  }
  # clean up any duplicate tip labels
  tree <- GetCleanTree(tree)
  return(tree)
}

# function to remove duplicated tip labels from a phylo or multiPhylo object
GetCleanTree <- function(tree){
  if(class(tree) == "phylo"){
    duptips <- which(duplicated(tree$tip.label))
    if (length(duptips) > 0)
      tree <- drop.tip(tree, duptips)
  }
  if(class(tree) == "multiPhylo"){
    for(i in 1:length(tree)){
      tr <- tree[[i]]
      duptips <- which(duplicated(tr$tip.label))
      if(length(duptips) > 0)
        tree[[i]] <- drop.tip(tr, duptips)
    }
  }
  return(tree)
}

# clean haploid values
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  for(i in 1:nrow(df)){
    x <- df$haploid[i]
    if(is.numeric(x)){
      if(!x == round(x)){
        df$haploid[i] <- stoch_round(x)
      }
    } else {
      nums <- as.numeric(strsplit(x, split="-", fixed=T)[[1]])
      if(length(nums) > 2){
        stop(paste("Looks like you have too many dashes! Check row:",i))
      }
      if(length(nums) == 2){
        df$haploid[i] <- sample(nums[1]:nums[2], 1)
      }
      nums <- as.numeric(strsplit(x, split=",", fixed=T)[[1]])
      if(length(nums)>1){
        df$haploid[i] <- sample(nums, 1)
      }
    }
  }
  df$haploid <- stoch_round(as.numeric(df$haploid))
  return(df)
}

# makes the tree ultrametric
AdjustTree <- function(tree){
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

###### End Functions ######

tree_dir  <- "../../data/trees/"
chrom_dir <- "../../data/chrome/"
# Get all files
tree.files <- list.files(tree_dir, full.names = FALSE)
karyo.files <- list.files(chrom_dir, full.names = FALSE)
clades <- tools::file_path_sans_ext(basename(tree.files))

i <- 1

for(i in 1:length(clades)){
  clade <- clades[i]
  print(paste("Working on", clade))
  # list of plant clades we run this because if you have one of these clades
  # we have to use slightly different functions later.
  plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "liliaceae", "solanaceae")
  # read chromosome data
  dat <- read.csv(paste0(chrom_dir, clade, ".csv"))
  dat <- CleanHaploid(dat)
  # read tree
  tree <- GetTree(clade)
  # Prune the phylogeny to only those species that you have chromosome data for. 
  if(class(tree) == "multiPhylo") tree <- tree[[1]]
  class(tree) <- "phylo"
  matched <- intersect(tree$tip.label, dat$species)
  tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
  
  # Prune chromosome data to only those that are present on the phylogeny.
  dat <- dat[dat$species %in% matched, ]
  dat <- dat[, c("species", "haploid")]
  dat <- dat[!duplicated(dat[,1]), ]
  # If trees are not ultrametric then we have to make them ultrametric
  if("phylo" %in% class(tree)){
    tree <- AdjustTree(tree)
  }
  vals <- dat$haploid
  names(vals) <- dat$species
  # 2. Setup the object
  obj <- contMap(tree, vals, plot = FALSE)
  obj <- setMap(obj, hcl.colors(n = 100, palette = "viridis"))
  
  # 1. Setup layout variables
  N <- length(obj$tree$tip.label)
  
  # 2. Plot the tree component manually
  #    ylim: Keeps the extra space at the bottom so the legend fits.
  #    ftype="off": Keeps the optimization crash fix.
  if(N<50) lwd <- 4
  if(N<100) lwd <- 3
  if(N<300) lwd <- 1
  if(N>1000) lwd <- .5
  png(paste(clade, "CM.png"), width = 5, height = 5, units = "in", res = 600)
  plot(obj$tree, 
       colors = obj$cols, 
       ylim = c(-0.15 * N, N), 
       ftype = "off", 
       outline = FALSE, 
       lwd = lwd)
  
  # 3. Add the legend with 'subtitle' suppressed
  #    We explicitly set subtitle="" to overwrite any default text below the bar.
  add.color.bar(leg = 0.5 * max(nodeHeights(obj$tree)), 
                cols = obj$cols, 
                title = "Haploid Count", 
                lims = obj$lims, 
                prompt = FALSE,
                x = 0, 
                y = -0.05 * N,
                subtitle = "")  
  dev.off()
  #i <- i + 1
}
