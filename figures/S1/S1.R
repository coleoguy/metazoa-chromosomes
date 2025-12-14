# visualize pull of the prior
# goal of this visualization is to determine the way that tree size 
# differences with exponential prior lead to biases in parameter estimates.

load(file = "../../results/pull.of.prior/simulation.analysis.RData")

library(ggplot2)
library(dplyr) 

n <- 100

# -----------------------------------------------------------------------------
# 1. DATA PROCESSING
# -----------------------------------------------------------------------------

# Add IDs to track replicates
hpd.full$id <- 1:nrow(hpd.full)
hpd.pruned$id <- 1:nrow(hpd.pruned)
sig.runs$id <- 1:nrow(sig.runs)

# Create a sorting key based on the 'Full Tree' lower bound
sort_order <- hpd.full[order(hpd.full$asc.low), "id"]

# Format "Full Tree" data
df_full <- data.frame(
  id = rep(hpd.full$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  tree_type = "Full Tree",
  low = c(hpd.full$asc.low, hpd.full$desc.low),
  high = c(hpd.full$asc.high, hpd.full$desc.high)
)

# Format "Pruned Tree" data
df_pruned <- data.frame(
  id = rep(hpd.pruned$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  tree_type = "Pruned Tree",
  low = c(hpd.pruned$asc.low, hpd.pruned$desc.low),
  high = c(hpd.pruned$asc.high, hpd.pruned$desc.high)
)

# Combine datasets
plot_data <- rbind(df_full, df_pruned)

# Format Significance data
sig_long <- data.frame(
  id = rep(sig.runs$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  is_significant = c(sig.runs$asc, sig.runs$desc)
)

# Merge and Sort
final_df <- merge(plot_data, sig_long, by = c("id", "parameter"))
final_df$id <- factor(final_df$id, levels = sort_order)

# -----------------------------------------------------------------------------
# 2. PLOTTING (Small Figure Optimized)
# -----------------------------------------------------------------------------

# Professional palette: Blue vs Red for clear contrast
cols <- c("Full Tree" = "#377EB8", "Pruned Tree" = "#E41A1C") 

ggplot(final_df, aes(y = id, color = tree_type)) +
  
  # GEOM LAYER
  geom_linerange(
    aes(xmin = low, xmax = high, 
        size = is_significant, 
        alpha = is_significant),
    position = position_identity()
  ) +
  
  # FACETING
  # EDITED: Removed "Parameter: " text
  facet_wrap(~parameter, scales = "free_x", 
             labeller = as_labeller(c(asc = "Ascending", 
                                      desc = "Descending"))) +
  
  # SCALES
  scale_color_manual(values = cols) +
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1.0), guide = "none") +
  # Compact line widths for small figure
  scale_size_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
  
  # THEME & LAYOUT
  theme_classic(base_size = 10) + 
  labs(
    x = "Parameter Estimate (95% HPD)",
    y = "Replicate (Sorted by Magnitude)",
    color = NULL 
  ) +
  
  theme(
    # FORCE WHITE BACKGROUND
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    
    # AXES
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    # REMOVE CLUTTER
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", margin = margin(t = 3)),
    
    # STRIPS (Facet Headers)
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10, hjust = 0),
    
    # LEGEND
    # EDITED: Moved to standard "top" position to prevent overlap
    legend.position = "top", 
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    
    # MARGINS
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )