# This script illustrates the importance of evaluating prior impact in cases
# where the goal is to compare across clades of different sizes.

library(chromePlus)
library(phytools)
library(ape)
library(TreeSim)
library(coda)
library(diversitree)

# simulate full trees
trees <- sim.bd.taxa(n = 500, numbsim = 100,
                     lambda = 3, mu = 1,
                     complete = FALSE)

# make trees unit length
dat <- res.full <- res.pruned <- list()
for(i in 1:100){
  # make trees unit length
  trees[[i]]$edge.length <- trees[[i]]$edge.length / max(branching.times(trees[[i]]))
  # simulate chromosome data
  dat[[i]] <- simChrom(tree = trees[[i]],pars=c(1,1, 0, 0, 20), limits = c(10, 30), model = "2010")
  #format for ChromePlus
  curdat <- data.frame(tips=names(dat[[i]]), haploid=dat[[i]])
  mat <- datatoMatrix(x = curdat, buffer = 1, hyper = FALSE)
  # make the likelihood function
  lik <- make.mkn(tree = trees[[i]], states = mat, k=ncol(mat), 
                  strict=FALSE,
                  control=list(method="ode", root=ROOT.OBS))
  conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                         polyploidy = FALSE, verbose=FALSE, constrain=list(drop.poly=T, drop.demi=T))
  print(argnames(conlik))
  # fit full data
  res.full[[i]] <- diversitree::mcmc(lik=conlik,
                        x.init=runif(length(argnames(conlik))),
                        prior=make.prior.exponential(2),
                        nsteps=500, w=1)
  # now pruned trees
  tree <- drop.tip(phy=trees[[i]], tip=sample(1:500, size=450, replace=F))
  curdat <- data.frame(tips=names(dat[[i]]), haploid=dat[[i]])
  curdat <- curdat[curdat$tips %in% tree$tip.label,]
  mat <- datatoMatrix(x = curdat, buffer = 1, hyper = FALSE)
  lik <- make.mkn(tree = tree, states = mat, k=ncol(mat), 
                  strict=FALSE,
                  control=list(method="ode", root=ROOT.OBS))
  conlik <- constrainMkn(data = mat, lik = lik, hyper=FALSE,
                         polyploidy = FALSE, verbose=FALSE, constrain=list(drop.poly=T, drop.demi=T))
  print(argnames(conlik))
  res.pruned[[i]] <- diversitree::mcmc(lik=conlik,
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
for(i in 1:100){
  hpd.full[i, 1:2] <- HPDinterval(as.mcmc(res.full[[i]]$asc1[10:100]))[1:2]
  hpd.full[i, 3:4] <- HPDinterval(as.mcmc(res.full[[i]]$desc1[10:100]))[1:2]
  hpd.pruned[i, 1:2] <- HPDinterval(as.mcmc(res.pruned[[i]]$asc1[10:100]))[1:2]
  hpd.pruned[i, 3:4] <- HPDinterval(as.mcmc(res.pruned[[i]]$desc1[10:100]))[1:2]
  if(pmax(hpd.full[i,1], hpd.pruned[i,1]) <= pmin(hpd.full[i,2], hpd.pruned[i,2])){
    sig.runs[i,1] <- F
  }else{
    sig.runs[i,1] <- T
  }
  if(pmax(hpd.full[i,3], hpd.pruned[i,3]) <= pmin(hpd.full[i,4], hpd.pruned[i,4])){
    sig.runs[i,2] <- F
  }else{
    sig.runs[i,2] <- T
  }
}
rm(list=ls()[-c(11,12,21,22,25)])



