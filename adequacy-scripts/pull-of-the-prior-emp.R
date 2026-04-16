# This script illustrates the importance of evaluating prior impact in cases
# where the goal is to compare across clades of different sizes.
library(ape)
library(chromePlus)
library(diversitree)
library(coda)
#library(phytools)
#library(TreeSim)
#

# Get Scarab Data
tree <- read.tree("../data/trees/scarabidae.new")[[1]]
as.phylo(tree) -> tree
chroms <- read.csv("../data/chrome/scarabidae.csv")[,3:4]
# prune to overlapping data
tree <- drop.tip(tree, tip=tree$tip.label[!tree$tip.label %in% chroms$species])
chroms <- chroms[chroms$species %in% tree$tip.label,]
# randomly select a haploid number for species that have multiple records
# 1. Split into a list of dataframes by species
species_list <- split(chroms, chroms$species)
# 2. Pick one random row from each (using sample)
#    We create a function that takes a dataframe 'df' and picks one row
sampled_list <- lapply(species_list, function(df) {
  df[sample(nrow(df), 1), ]
})
# 3. Combine back into one dataframe
chroms <- do.call(rbind, sampled_list)

# make trees unit length
tree$edge.length <- tree$edge.length / max(branching.times(tree))

mat <- datatoMatrix(x = chroms, buffer = 1, hyper = FALSE)
# make the likelihood function
lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                strict=FALSE,
                control=list(method="ode", root=ROOT.OBS))
conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                       polyploidy = FALSE, verbose=FALSE, constrain=list(drop.poly=T, drop.demi=T))
print(argnames(conlik))
# fit full data
res.full <- diversitree::mcmc(lik=conlik,
                      x.init=runif(length(argnames(conlik))),
                      prior=make.prior.exponential(2),
                      nsteps=500, w=1)
res.pruned <- list()

for(i in 1:100){
  curtree <- drop.tip(tree, tip=sample(1:length(tree$tip.label), size=round(.5*length(tree$tip.label)), replace=F))

  curdat <- chroms[chroms$species %in% curtree$tip.label,]
  mat <- datatoMatrix(x = curdat, buffer = 1, hyper = FALSE)
  # make the likelihood function
  lik <- make.mkn(tree = curtree, states = mat, k=ncol(mat), 
                  strict=FALSE,
                  control=list(method="ode", root=ROOT.OBS))
  conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                         polyploidy = FALSE, verbose=FALSE, constrain=list(drop.poly=T, drop.demi=T))
  print(argnames(conlik))
  # fit full data
  res.pruned[[i]] <- mcmc(lik=conlik,
                        x.init=runif(length(argnames(conlik))),
                        prior=make.prior.exponential(2),
                        nsteps=500, w=1)
}


### move to a plotting file ###
hpd.full <- hpd.pruned <- data.frame(asc.low = NA,
                                     asc.high = NA,
                                     desc.low = NA,
                                     desc.high =NA)
sig.runs <- data.frame(asc = NA,
                       desc = NA)
hpd.full[1, 1:2] <- HPDinterval(as.mcmc(res.full$asc1[10:100]))[1:2]
hpd.full[1, 3:4] <- HPDinterval(as.mcmc(res.full$desc1[10:100]))[1:2]



for(i in 1:100){
  hpd.pruned[i, 1:2] <- HPDinterval(as.mcmc(res.pruned[[i]]$asc1[10:100]))[1:2]
  hpd.pruned[i, 3:4] <- HPDinterval(as.mcmc(res.pruned[[i]]$desc1[10:100]))[1:2]
  if(pmax(hpd.full[1,1], hpd.pruned[i,1]) <= pmin(hpd.full[1,2], hpd.pruned[i,2])){
    sig.runs[i,1] <- F
  }else{
    sig.runs[i,1] <- T
  }
  if(pmax(hpd.full[1,3], hpd.pruned[i,3]) <= pmin(hpd.full[1,4], hpd.pruned[i,4])){
    sig.runs[i,2] <- F
  }else{
    sig.runs[i,2] <- T
  }
}
rm(list=ls()[-c(5,6,10,11,13)])



