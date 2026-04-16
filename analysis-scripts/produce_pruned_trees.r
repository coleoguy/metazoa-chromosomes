### Megan Copeland

## Load libraries
library(ape)           
library(diversitree)   
library(chromePlus)
library(phangorn)

###### Functions ######
# function to read and return a tree object for a given clade name
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0("../data/trees/", x, ".new")  # possible Newick file path
  nex_path <- paste0("../data/trees/", x, ".nex")  # possible Nexus file path
  
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
  if(class(tree) == "phylo"){                   # if single tree
    duptips <- which(duplicated(tree$tip.label))   # find duplicate tips
    if (length(duptips) > 0)
      tree <- drop.tip(tree, duptips)              # remove them
  }
  if(class(tree) == "multiPhylo"){       # if multiple trees
    for(i in 1:length(tree)){
      tr <- tree[[i]]
      duptips <- which(duplicated(tr$tip.label))
      if(length(duptips) > 0)
        tree[[i]] <- drop.tip(tr, duptips)
    }
  }
  return(tree)
}
# function to clean haploid column (to handle ranges, comma seperated values, multiple entries)
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  for(i in 1:nrow(df)){
    x <- df$haploid[i]
    if(is.numeric(x)){
      if(!x == round(x)){
        df$haploid[i] <- stoch_round(x)
      }
    }else{
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

AdjustTree <- function(tree){
  # ensure the tree is ultrametric
  if (!is.ultrametric(tree)) {
    if (clade %in% plants) {
      tree <- nnls.tree(cophenetic(tree), tree, rooted = T)
    } else {
      tree <- chronos(tree)
    }
  }
  # scale tree to unit length
  tree$edge.length <- tree$edge.length / max(branching.times(tree))
  # replace any small negative branch lengths with tiny positive value
  tree$edge.length[tree$edge.length < 0] <- 1e-9
  # ensure tree is fully bifurcating 
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
  }
  return(tree)
}


###### End Functions ######

## Main Loop

# get all plant clades that are using the angiosperm phylo
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae")

# list all chromosome data files
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)

# initialize an empty list to store MCMC results
results <- list()

# iterate through each chrome data file
for (i in 1:length(file_list)) {
  f <- file_list[[i]]
  clade <- tools::file_path_sans_ext(basename(f))  # extract clade name from filename
  print(clade)
  
  # read chromosome data for this clade
  dat <- read.csv(f)
  dat <- CleanHaploid(dat)
  
  # read the corresponding tree
  tree <- GetTree(clade)
  
  # match tree and data by species names
  if(class(tree) == "phylo"){
    matched <- intersect(tree$tip.label, dat$species)
    tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
  }
  if(class(tree) == "multiPhylo"){
    result <- list()
    for(i in 1:length(tree)){
      tr <- tree[[i]]
      matched <- intersect(tr$tip.label, dat$species)
      result[[i]] <- keep.tip(tr, matched)
    }
    tree <- result
    class(tree) <- "multiPhylo"
  }
  write.nexus(tree, file=paste("../data/pruned-trees/", clade, "pruned.nex", sep=""))
}
