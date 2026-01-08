### Megan Copeland and Heath Blackmon

## Load libraries
library(parallel)
library(doParallel)
library(foreach)

library(ape)
library(diversitree)
library(chromePlus)
library(phangorn)

# 1. Detect cores
num_cores <- 20

# 2. Setup the cluster
my_cluster <- makeCluster(num_cores)

# 3. Register the cluster
registerDoParallel(my_cluster)

# list all chromosome data files
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = TRUE)

# Parallel Loop
parallel_results <- foreach(i = 1:length(file_list), .combine = 'list') %dopar% {
  
  ## Load libraries inside worker
  library(ape)
  library(diversitree)
  library(chromePlus)
  library(phangorn)
  
  ###### Functions (MATCH YOUR "GOOD" SCRIPT EXACTLY) ######
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
  
  ###### Settings ######
  plants <- c("asteraceae", "fabaceae", "brassicaceae",
              "orchidaceae", "lilaceae", "solanaceae")
  
  burnin <- 500
  target_params <- c("asc1", "desc1", "pol1", "dem1")
  
  model_names <- c("full", "no_pol", "no_dem", "no_poly_no_demi")
  ###### Settings ######
  
  ###### Identify clade + paths ######
  f <- file_list[[i]]
  clade <- tools::file_path_sans_ext(basename(f))
  
  ex.result <- read.csv("../results/AIC.csv")
  
  if(!clade %in% ex.result$clade){
    
    mcmc_path <- paste0("../results/exponential.prior - full model/mentor_results/",
                        clade, ".p.res.csv")
    
    ###### Get start values from MCMC ######
    dat_mcmc <- read.csv(mcmc_path, check.names = FALSE)
    keep <- dat_mcmc$i > burnin
    
    mcmc_means <- rep(NA_real_, length(target_params))
    names(mcmc_means) <- target_params
    
    for(p in target_params){
      if(p %in% names(dat_mcmc)){
        vals <- dat_mcmc[keep, p]
        vals <- vals[is.finite(vals) & vals > 0]
        if(length(vals) > 0)
          mcmc_means[p] <- mean(vals)
      }
    }
    
    ###### Read chromosome data (CLEAN FIRST, LIKE YOUR GOOD SCRIPT) ######
    dat_full <- read.csv(f, stringsAsFactors = FALSE)
    dat_full <- CleanHaploid(dat_full)
    
    ###### Read tree (AFTER CLEANING, LIKE YOUR GOOD SCRIPT) ######
    tree <- GetTree(clade)
    if(is.null(tree)){
      return(clade)
    }
    
    ###### Prune tree to taxa w/ chromosome data (MATCH YOUR GOOD SCRIPT STYLE) ######
    if(class(tree) == "phylo"){
      matched <- intersect(tree$tip.label, dat_full$species)
      tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
    }
    if(class(tree) == "multiPhylo"){
      result <- list()
      for(j in 1:length(tree)){
        tr <- tree[[j]]
        matched <- intersect(tr$tip.label, dat_full$species)
        result[[j]] <- keep.tip(tr, matched)
      }
      tree <- result
      class(tree) <- "multiPhylo"
    }
    
    dat <- dat_full[dat_full$species %in% matched, ]
    dat <- dat[, c("species", "haploid")]
    
    mat <- datatoMatrix(x = dat, buffer = 1, hyper = FALSE)
    
    if("phylo" %in% class(tree)){
      tree <- AdjustTree(tree)
    }
    if("multiPhylo" %in% class(tree)){
      for(j in 1:length(tree)){
        tree[[j]] <- AdjustTree(tree[[j]])
        print(paste("working on tree", j))
      }
    }
    
    ###### models ######
    models <- list(
      full = function(mat, lik) {
        constrainMkn(data = mat, lik = lik, hyper=FALSE,
                     polyploidy = FALSE, verbose=FALSE)
      },
      no_pol = function(mat, lik) {
        constrainMkn(data = mat, lik = lik, hyper=FALSE,
                     polyploidy = FALSE, verbose=FALSE,
                     constrain = list(drop.poly = TRUE))
      },
      no_dem = function(mat, lik) {
        constrainMkn(data = mat, lik = lik, hyper=FALSE,
                     polyploidy = FALSE, verbose=FALSE,
                     constrain = list(drop.demi = TRUE))
      },
      no_poly_no_demi = function(mat, lik) {
        constrainMkn(data = mat, lik = lik, hyper=FALSE,
                     polyploidy = FALSE, verbose=FALSE,
                     constrain = list(drop.poly = TRUE, drop.demi = TRUE))
      }
    )
    
    ###### Run models (per tree) ######
    clade_aic_results <- list()
    
    if("phylo" %in% class(tree)){
      
      lik <- make.mkn(tree = tree, states = mat, k=ncol(mat),
                      strict=FALSE,
                      control=list(method="ode", root=ROOT.OBS))
      
      tree_row <- c(full=NA, no_pol=NA, no_dem=NA, no_poly_no_demi=NA)
      
      for(m in model_names){
        conlik_m <- models[[m]](mat, lik)
        par_m <- argnames(conlik_m)
        x0 <- mcmc_means[par_m]
        
        if(all(is.finite(x0)) && all(x0 > 0)){
          fit <- tryCatch(
            find.mle(conlik_m, x.init = x0, method = "subplex", lower = 0, upper = 30),
            error = function(e) e
          )
          if(!inherits(fit, "error"))
            tree_row[m] <- AIC(fit)
        }
      }
      
      clade_aic_results[[1]] <- tree_row
    }
    
    if("multiPhylo" %in% class(tree)){
      
      for(j in 1:length(tree)){
        
        tr <- tree[[j]]
        
        # NOTE: we keep the same mat for all trees because dat was pruned once
        lik <- make.mkn(tree = tr, states = mat, k=ncol(mat),
                        strict=FALSE,
                        control=list(method="ode", root=ROOT.OBS))
        
        tree_row <- c(full=NA, no_pol=NA, no_dem=NA, no_poly_no_demi=NA)
        
        for(m in model_names){
          
          conlik_m <- models[[m]](mat, lik)
          par_m <- argnames(conlik_m)
          x0 <- mcmc_means[par_m]
          
          if(all(is.finite(x0)) && all(x0 > 0)){
            fit <- tryCatch(
              find.mle(conlik_m, x.init = x0, method = "subplex", lower = 0, upper = 30),
              error = function(e) e
            )
            if(!inherits(fit, "error"))
              tree_row[m] <- AIC(fit)
          }
        }
        
        clade_aic_results[[j]] <- tree_row
      }
    }
    
    ###### Aggregate + write results ######
    res_mat <- do.call(rbind, clade_aic_results)
    mean_aics <- colMeans(res_mat, na.rm = TRUE)
    out_df <- data.frame(t(mean_aics))
    rownames(out_df) <- clade
    write.csv(out_df, file=paste0("../results/model.testing/", clade, "_AIC.csv"))
  }
    
   return(clade)
}

# Always stop cluster
stopCluster(my_cluster)


#write all model testing results out to file
#files <- list.files("../results/model.testing (some missing results)", full.names = TRUE)
#dat <- do.call(rbind,  lapply(files, read.csv))
#colnames(dat) <- c("clade", colnames(dat[2:5]))
#write.csv(dat, "../results/AIC.csv", row.names = F)
