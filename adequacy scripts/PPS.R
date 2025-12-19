## Load libraries
library(ape)
library(diversitree)
library(chromePlus)
library(phangorn)

###### Functions ######

# Reuse your existing tree loading logic
GetTree <- function(x){
  tree <- NULL
  new_path <- paste0("../data/trees/", x, ".new")
  nex_path <- paste0("../data/trees/", x, ".nex")
  
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

GetCleanTree <- function(tree){
  if(class(tree) == "phylo"){
    duptips <- which(duplicated(tree$tip.label))
    if (length(duptips) > 0) tree <- drop.tip(tree, duptips)
  }
  if(class(tree) == "multiPhylo"){
    for(i in 1:length(tree)){
      tr <- tree[[i]]
      duptips <- which(duplicated(tr$tip.label))
      if(length(duptips) > 0) tree[[i]] <- drop.tip(tr, duptips)
    }
  }
  return(tree)
}

# Reuse your cleaning function for the Empirical Data
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  for(i in 1:nrow(df)){
    x <- df$haploid[i]
    if(is.numeric(x)){
      if(!x == round(x)) df$haploid[i] <- stoch_round(x)
    } else {
      nums <- as.numeric(strsplit(x, split="-", fixed=T)[[1]])
      if(length(nums) == 2) df$haploid[i] <- sample(nums[1]:nums[2], 1)
      nums <- as.numeric(strsplit(x, split=",", fixed=T)[[1]])
      if(length(nums)>1) df$haploid[i] <- sample(nums, 1)
    }
  }
  df$haploid <- stoch_round(as.numeric(df$haploid))
  return(df)
}

# Helper to calculate Variance and Shannon's Entropy
GetStats <- function(chrom_vals) {
  # remove NAs if any
  chrom_vals <- na.omit(chrom_vals)
  
  # Variance
  val_var <- var(chrom_vals)
  
  # Shannon's Entropy
  # H = -sum(p_i * ln(p_i))
  counts <- table(chrom_vals)
  probs <- counts / sum(counts)
  val_ent <- -sum(probs * log(probs))
  
  return(c(variance = val_var, entropy = val_ent))
}

###### End Functions ######
n_sims <- 1000              # Number of simulations
burnin_prop <- 0.25         # Proportion of MCMC to discard as burnin
file_list <- list.files(path = "../data/chrome", pattern = "\\.csv$", full.names = F)
clades <- unlist(strsplit(file_list, ".", fixed=T))[seq(from=1, by=2, length.out= 56)]
plants <- c("asteraceae", "fabaceae", "brassicaceae", "orchidaceae", "lilaceae", "solanaceae")

for(m in 48:length(clades)){
  
  
  ## ----- CONFIGURATION -----
  clade <- clades[m]      # <-- choose your clade
  # Paths
  data_path <- paste0("../data/chrome/", clade, ".csv")
  mcmc_path <- paste0("../results/uniform prior/",clade, ".p.res.csv") # Assuming MCMC output is in current dir
  tree_file <- clade
  
  ## ----- LOAD & PREP DATA -----
  
  # 1. Load Tree
  tree <- GetTree(tree_file)
  if(is.null(tree)) stop("Tree file not found.")
  
  # 2. Load Empirical Data
  dat <- read.csv(data_path)
  dat <- CleanHaploid(dat)
  
  # 3. Prune Tree/Data to match (Critical for accurate comparisons)
  if(class(tree) == "phylo"){
    matched <- intersect(tree$tip.label, dat$species)
    tree <- drop.tip(tree, setdiff(tree$tip.label, matched))
    # makes the tree ultrametric
    tree <- AdjustTree(tree)
    # Standardize edge lengths for chromePlus
    tree$edge.length <- tree$edge.length / max(branching.times(tree))
  }
  if("multiPhylo" %in% class(tree)){
    # If multiPhylo, we just prune the first one to get the species list
    # We assume all trees in the file have the same tip set
    matched <- intersect(tree[[1]]$tip.label, dat$species)
    foo <- list()
    for(i in 1:length(tree)){
      foo[[i]] <- keep.tip(tree[[i]], matched)
      foo[[i]] <- AdjustTree(foo[[i]])
      foo[[i]]$edge.length <- foo[[i]]$edge.length / max(branching.times(foo[[i]]))
    }
    tree <- foo
    class(tree) <- "multiPhylo"
  }
  
  # Filter data to matched species
  dat <- dat[dat$species %in% matched, ]
  
  # 4. Load MCMC Results
  mcmc <- read.csv(mcmc_path)
  
  # Remove Burnin
  burnin_rows <- floor(nrow(mcmc) * burnin_prop)
  mcmc_post <- mcmc[(burnin_rows + 1):nrow(mcmc), ]
  
  # Sample 1000 rows from the posterior for simulation
  # If we have fewer than 1000 post-burnin samples, we sample with replacement
  sample_indices <- sample(1:nrow(mcmc_post), n_sims, replace = (nrow(mcmc_post) < n_sims))
  mcmc_sample <- mcmc_post[sample_indices, ]
  
  ## ----- CALCULATE EMPIRICAL STATS -----
  emp_stats <- GetStats(dat$haploid)
  print(paste("Empirical Variance:", round(emp_stats['variance'], 3)))
  print(paste("Empirical Entropy:", round(emp_stats['entropy'], 3)))
  
  ## ----- RUN SIMULATIONS -----
  print("Starting Simulations...")
  
  # Initialize storage for results (1000 rows, 2 columns)
  sim_results <- matrix(NA, nrow = n_sims, ncol = 2)
  colnames(sim_results) <- c("variance", "entropy")
  
  # We define the range for simulation (limits). 
  # Good practice: 1 to max_observed + buffer (e.g., 20)
  max_chrom <- max(dat$haploid) + 5
  limits <- c(1, max_chrom)
  
  for(i in 1:n_sims){
    
    # 1. Get parameters for this simulation
    pars_row <- mcmc_sample[i,]
    
    # Extract rates. 
    # Note: chromePlus::simChrom expects a vector 'pars'. 
    # We assume the MCMC columns match the necessary names (asc1, desc1, etc.)
    pars_vector <- c()
    partypes <- c("asc1","desc1","dem1","pol1")
    # have <- unlist(pars_row[ , names(pars_row) %in% c("asc1","desc1","dem1","pol1")])
    pars_vector[1] <- pars_row[, names(pars_row) == "asc1"]
    pars_vector[2] <- pars_row[, names(pars_row) == "desc1"]
    if("dem1" %in% names(pars_row)){
      pars_vector[3] <- pars_row[, names(pars_row) == "dem1"]
    }else{
      pars_vector[3] <- 0
    }
    if("pol1" %in% names(pars_row)){
      pars_vector[4] <- pars_row[, names(pars_row) == "pol1"]
    }else{
      pars_vector[4] <- 0
    }
    pars_vector[5] <- round(mean(dat$haploid))
    
    # 2. Pick the correct tree
    # If single tree, use 'tree'. If multiPhylo, check if MCMC has 'tree' column.
    current_tree <- NULL
    if("phylo" %in% class(tree)){
      current_tree <- tree
    } else if ("multiPhylo" %in% class(tree)){
        current_tree <- tree[[sample(1:length(tree), 1)]]
      }
    
    
    # 3. Simulate
    # If your model assumes the root is the observed mean or fixed, simChrom usually samples root freq.
    # We use the parameters to simulate.
    tryCatch({
      
      
      
      sim_data <- simChrom(tree = current_tree, 
                           pars = pars_vector, 
                           limits = limits, 
                           model = "2010") # "2010" is standard standard Mk-style
      
      # 4. Calculate Stats on Simulation
      # simChrom returns a matrix where column 1 is tip names, col 2 is chrom number
      # (Check structure of sim_data, sometimes it returns just the vector depending on version)
      
      # If sim_data is a matrix/df:
      sim_vals <- sim_data
      
      sim_results[i, ] <- GetStats(sim_vals)
      
    }, error = function(e){
      message(paste("Simulation", i, "failed:", e$message))
      # Leave as NA
    })
    
    if(i %% 100 == 0) print(paste("Finished simulation", i, "of", n_sims))
  }
  
  ## ----- WRITE OUTPUT -----
  
  # Combine Empirical (row 1) and Simulated (rows 2-1001)
  final_output <- rbind(emp_stats, sim_results)
  row.names(final_output) <- c("empirical", paste0("sim_", 1:n_sims))
  
  write.csv(final_output, paste0("../results/PPS/",clade, ".PPS.csv"))
}
print("Done! Results written to CSV.")