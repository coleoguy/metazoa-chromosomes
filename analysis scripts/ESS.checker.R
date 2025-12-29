library(coda)

# Get the list of files
file.list <- list.files(path = "../results/exponential.prior/mentor_results", 
                        full.names = TRUE)

# Initialize an empty data frame to store results
# We include a 'file' column so you know which row belongs to which result
ess_table <- data.frame(
  file  = character(),
  asc1  = numeric(),
  desc1 = numeric(),
  demi  = numeric(),
  poly  = numeric(),
  stringsAsFactors = FALSE
)

for (f in file.list) {
  
  # 1. Read the file
  # We use tryCatch to skip files that might be empty or corrupted
  chain <- tryCatch(read.csv(f), error = function(e) NULL)
  
  if (!is.null(chain) && nrow(chain) > 100) {
    
    # 2. Remove 100 generations of burnin
    chain <- chain[-c(1:100), ]
    
    # 3. Define the parameters we want to measure
    params <- c("asc1", "desc1", "dem1", "pol1")
    
    # 4. Calculate ESS for each parameter if it exists
    # We create a named vector of NAs first
    row_values <- setNames(rep(NA, length(params)), params)
    
    for (p in params) {
      if (p %in% colnames(chain)) {
        # Check if the column has variance (ESS fails on constant values)
        if (var(chain[[p]]) > 0) {
          row_values[p] <- effectiveSize(as.mcmc(chain[[p]]))
        } else {
          row_values[p] <- 0 # or NA, depending on preference for fixed parameters
        }
      }
    }
    
    # 5. Add to the main table
    ess_table <- rbind(ess_table, data.frame(
      file  = basename(f),
      asc1  = row_values["asc1"],
      desc1 = row_values["desc1"],
      demi  = row_values["dem1"],
      poly  = row_values["pol1"]
    ))
  }
}

# save final table
write.csv(ess_table, file="../results/exponential.prior/ess_summary.csv", row.names = FALSE)
