## Load libraries
library(ape)           
library(diversitree)   
library(chromePlus)
library(phangorn)

###### Functions ######
# function to read and return a tree object for a given clade name
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0("../data/trees/", x, ".new")  # <-- change path as needed
  nex_path <- paste0("../data/trees/", x, ".nex")  # <-- change path as needed

  
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

## ----- SINGLE CLADE RUN -----
clade <- "solanaceae"   # <-- choose your clade

f <- paste0("../data/chrome/", clade, ".csv") # <-- change path as needed

print(clade)

# list of plant clades
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae")

# initialize results list
results <- list()

# read chromosome data
dat <- read.csv(f)
dat <- CleanHaploid(dat)

# read tree
tree <- GetTree(clade)

# match tree + data
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

# subset data
dat <- dat[dat$species %in% matched, ]
dat <- dat[, c("species", "haploid")]

# convert to diversitree format
mat <- datatoMatrix(x = dat, buffer = 1, hyper = FALSE)

# adjust tree(s)
if("phylo" %in% class(tree)){
  tree <- AdjustTree(tree)
}

if("multiPhylo" %in% class(tree)){
  for(j in 1:length(tree)){
    tree[[j]] <- AdjustTree(tree[[j]])
    print(paste("working on tree", j))
  }
}

# run model(s)
if("phylo" %in% class(tree)){
  lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                  strict=FALSE,
                  control=list(method="ode", root=ROOT.OBS))
  
  conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                         polyploidy = FALSE, verbose=FALSE)
  
  print(argnames(conlik))
  
  res <- mcmc(lik=conlik,
              x.init=runif(length(argnames(conlik))),
              prior=make.prior.exponential(2),
              nsteps=1000, w=1)
}

if("multiPhylo" %in% class(tree)){
  res <- list()
  for(j in 1:100){
    lik <- make.mkn(tree = tree[[j]], states = mat, k=ncol(mat),
                    strict=FALSE,
                    control=list(method="ode", root=ROOT.OBS))
    
    conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                           polyploidy = FALSE, verbose=FALSE)
    
    print(argnames(conlik))
    
    res[[j]] <- mcmc(lik=conlik,
                     x.init=runif(length(argnames(conlik))),
                     prior=make.prior.exponential(2),
                     nsteps=1000, w=1)
  }
}

## ----- WRITE OUT RESULTS -----

# single tree results
if("phylo" %in% class(tree)){

  # write csv
  write.csv(as.data.frame(res),
            file = paste0(clade, "_mcmc.csv"),
            row.names = FALSE)
}

# multiple tree results
if("multiPhylo" %in% class(tree)){

  # loop in the same style as earlier code
  for(j in 1:length(res)){
    this <- res[[j]]
    
    # skip if NULL (failed run)
    if(!is.null(this)){
      write.csv(as.data.frame(this),
                file = paste0(clade, "_tree", j, ".csv"),
                row.names = FALSE)
    }
  }
}
