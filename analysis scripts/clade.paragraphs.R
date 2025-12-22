library(ape)
library(coda)

# --- 1. Setup ---
higher_df  <- read.csv("../figures/figure 1/higher_class.csv")
mcmc_files <- list.files("../results/exponential.prior/mentor_results", full.names = TRUE)
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
GetPolytomyStats <- function(phy){
  # 1. Basic Setup
  n_tips <- Ntip(phy)
  parents <- phy$edge[, 1]
  children <- phy$edge[, 2]
  # 2. Identify Polytomous Nodes (nodes with degree > 2)
  # table() counts the occurrence of each parent node
  node_counts <- table(parents)
  poly_nodes <- as.numeric(names(node_counts)[node_counts > 2])
  # 3. Identify Tips in Polytomies
  # Logic: The child must be a tip (index <= n_tips) AND
  #        The parent must be in the list of polytomous nodes
  tips_in_polytomy <- sum((children <= n_tips) & (parents %in% poly_nodes))
  # 4. Calculate Percentage
  percent <- (tips_in_polytomy / n_tips) * 100
  # 5. Return Vector
  return(c(count = tips_in_polytomy, percent = percent))
}

extdat <- data.frame(
  clade = NA,
  h.taxonomy = NA,
  chromes = NA,
  phytips = NA,
  overlap = NA,
  rootage = NA,
  taxtips = NA,
  unresolved.perc = NA)

# Higher DF will determine the order of the paragraphs
for(i in 1:length(tree$tip.label)){
  curclade <- tree$tip.label[i]
  chromes  <- read.csv(paste0("../data/chrome/", curclade, ".csv"))
  curtree <- GetTree(curclade)
  if(class(curtree) == "multiPhylo"){
    curtree <- curtree[[1]]
  }
  extdat[i,1] <- curclade
  extdat[i,2] <- higher_df$Higher.Classification[higher_df$Clade == tolower(curclade)]
  extdat[i,3] <- nrow(chromes)
  extdat[i,4] <- length(curtree$tip.label)
  keep <- curtree$tip.label[curtree$tip.label %in% chromes$species]
  p.tree <- keep.tip(curtree, keep)
  
  extdat[i,5] <- length(p.tree$tip.label)
  extdat[i,6] <- ages$Age[ages$Clade == curclade]
  extdat[i,7] <- !is.binary(p.tree)
  if(extdat[i,7]){
    extdat[i,8] <- GetPolytomyStats(p.tree)[2]
  }
}
write.csv(extdat, file="../results/data.stats.csv", row.names = F)



