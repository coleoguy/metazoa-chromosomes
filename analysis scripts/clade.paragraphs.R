library(ape)
library(coda)

# --- 1. Setup ---
higher_df  <- read.csv("../figures/figure 1/higher_class.csv")
mcmc_files <- list.files("../results/exponential.prior/mentor_results/", full.names = TRUE)
ages       <- read.csv("../data/clade.ages.csv")
tree       <- read.tree("../figures/figure 1/cladetree.new")
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

extdat <- data.frame(
  clade = NA,
  h.taxonomy = NA,
  chromes = NA,
  phytips = NA,
  overlap = NA,
  rootage = NA,
  taxtips = NA)

# Higher DF will determine the order of the paragraphs
for(i in length(tree$tip.label):1){
  curclade <- tree$tip.label[i]
  chromes  <- read.csv(paste0("../data/chrome/", curclade, ".csv"))
  curtree <- GetTree(curclade)
  extdat$clade[i] <- curclade
  extdat$h.taxonomy[i] <- higher_df$Higher.Classification[higher_df$Clade == curclade]
  extdat$chromes [i] <- nrow(chromes)
  extdat$phytips [i] <- tree$tip.label[i]
  extdat$overlap [i] <- tree$tip.label[i]
  extdat$rootage [i] <- tree$tip.label[i]
  extdat$taxtips [i] <- tree$tip.label[i]
  extdat$range [i] <- tree$tip.label[i]

}




