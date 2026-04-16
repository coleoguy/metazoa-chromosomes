# -----------------------------------------------------------------------------
# 1. SETUP: Define your model folders here
# Format: "Model_Name" = "Path/To/Folder"
# -----------------------------------------------------------------------------
model_folders <- c(
  "With_Poly"    = "../results/polyploidy_model",
  "Without_Poly" = "../results/simple_model"
)

# Initialize the final merged dataframe as NULL
final_results <- NULL

# -----------------------------------------------------------------------------
# 2. LOOP: Process each model folder
# -----------------------------------------------------------------------------
for (model_name in names(model_folders)) {
  
  dir_path <- model_folders[[model_name]]
  file_list <- list.files(dir_path, full.names = TRUE)
  
  # Vectors to temporarily store data for this specific model
  clade_names <- character()
  median_rates <- numeric()
  
  # Inner loop: Process each file in the current folder
  for (f in file_list) {
    
    # Read file safely
    # We use tryCatch to skip corrupt/empty files without crashing
    dat <- tryCatch(read.csv(f), error = function(e) NULL)
    
    # Check if data exists, has enough rows for burnin, and has required columns
    if (!is.null(dat) && nrow(dat) > 100 && all(c("asc1", "desc1") %in% colnames(dat))) {
      
      # Remove first 100 rows (Burnin)
      dat <- dat[-c(1:100), ]
      
      # Calculate sum of rows, then median
      # Row-wise sum in base R is just vec1 + vec2
      current_median <- median(dat$asc1 + dat$desc1)
      
      # Clean clade name (remove extension)
      # Adjust pattern if your files are .rds or something else
      current_clade <- sub("\\.[^.]+$", "", basename(f))
      
      # Append to vectors
      clade_names <- c(clade_names, current_clade)
      median_rates <- c(median_rates, current_median)
    }
  }
  
  # Create a dataframe for this model
  model_df <- data.frame(
    clade = clade_names,
    rate = median_rates,
    stringsAsFactors = FALSE
  )
  
  # Calculate Rank (1 = Highest Rate)
  # rank(-x) gives rank 1 to the largest number
  model_df$rank <- rank(-model_df$rate, ties.method = "min")
  
  # Rename columns to include the model name (e.g., rate.With_Poly, rank.With_Poly)
  colnames(model_df)[2:3] <- paste(c("rate", "rank"), model_name, sep = ".")
  
  # Merge into the final dataframe
  if (is.null(final_results)) {
    final_results <- model_df
  } else {
    # Merge by "clade", keeping all rows (all = TRUE behaves like full_join)
    final_results <- merge(final_results, model_df, by = "clade", all = TRUE)
  }
}

# -----------------------------------------------------------------------------
# 3. OUTPUT: View or Save
# -----------------------------------------------------------------------------
# Sort by the first model's rank for readability
if (ncol(final_results) >= 3) {
  first_rank_col <- grep("rank", colnames(final_results))[1]
  final_results <- final_results[order(final_results[, first_rank_col]), ]
}

# Print the top of the table
head(final_results)

# Optional: Write to CSV
# write.csv(final_results, "model_comparison_ranks.csv", row.names = FALSE)