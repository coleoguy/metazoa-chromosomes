# visualize pull of the prior
# Comparison: Single Full Tree Run vs. 100 Pruned Tree Replicates

load(file = "../../results/pull.of.prior/empirical.analysis.RData")

library(ggplot2)
library(dplyr) 

n <- 100

# -----------------------------------------------------------------------------
# 1. DATA WRANGLING (Base R)
# -----------------------------------------------------------------------------

# NOTE: Assuming hpd.full has 1 row and hpd.pruned has 100 rows.
# We replicate the single hpd.full row 100 times to create the visual "band".
hpd.full.expanded <- hpd.full[rep(1, n), ]

# Recalculate Significance (Overlap check)
# Significant = No Overlap between the specific Pruned run and the single Full run
# Logic: (Low_Pruned > High_Full) OR (High_Pruned < Low_Full)
sig.runs <- data.frame(
  asc = (hpd.pruned$asc.low > hpd.full.expanded$asc.high) | 
    (hpd.pruned$asc.high < hpd.full.expanded$asc.low),
  desc = (hpd.pruned$desc.low > hpd.full.expanded$desc.high) | 
    (hpd.pruned$desc.high < hpd.full.expanded$desc.low)
)

# Add IDs
hpd.full.expanded$id <- 1:n
hpd.pruned$id <- 1:n
sig.runs$id <- 1:n

# Create Sorting Key
# We MUST sort by 'Pruned Tree' now, because 'Full Tree' is constant.
# Sorting by pruned lower bound creates the caterpillar S-curve.
sort_order <- hpd.pruned[order(hpd.pruned$asc.low), "id"]

# --- Format "Full Tree" (Reference) ---
df_full <- data.frame(
  id = rep(hpd.full.expanded$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  tree_type = "Full Tree",
  low = c(hpd.full.expanded$asc.low, hpd.full.expanded$desc.low),
  high = c(hpd.full.expanded$asc.high, hpd.full.expanded$desc.high)
)

# --- Format "Pruned Tree" (Replicates) ---
df_pruned <- data.frame(
  id = rep(hpd.pruned$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  tree_type = "Pruned Tree",
  low = c(hpd.pruned$asc.low, hpd.pruned$desc.low),
  high = c(hpd.pruned$asc.high, hpd.pruned$desc.high)
)

# --- Combine ---
plot_data <- rbind(df_full, df_pruned)

# --- Format Significance ---
sig_long <- data.frame(
  id = rep(sig.runs$id, 2),
  parameter = c(rep("asc", n), rep("desc", n)),
  is_significant = c(sig.runs$asc, sig.runs$desc)
)

# --- Merge and Apply Sort Order ---
final_df <- merge(plot_data, sig_long, by = c("id", "parameter"))
final_df$id <- factor(final_df$id, levels = sort_order)

# -----------------------------------------------------------------------------
# 2. PLOTTING (Small Figure Optimized)
# -----------------------------------------------------------------------------

cols <- c("Full Tree" = "#377EB8", "Pruned Tree" = "#E41A1C") 

ggplot(final_df, aes(y = id, color = tree_type)) +
  
  # GEOM LAYER
  # This will visually create a vertical blue "band" (the constant Full Tree)
  # intersected by the varying red lines (Pruned Tree).
  geom_linerange(
    aes(xmin = low, xmax = high, 
        size = is_significant, 
        alpha = is_significant),
    position = position_identity()
  ) +
  
  # FACETING
  facet_wrap(~parameter, scales = "free_x", 
             labeller = as_labeller(c(asc = "Ascending", 
                                      desc = "Descending"))) +
  
  # SCALES
  scale_color_manual(values = cols) +
  # Alpha: Overlap = Transparent (mixes colors), Significant = Opaque
  scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1.0), guide = "none") +
  # Size: Thin lines for dense plot, slightly thicker for significant hits
  scale_size_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
  
  # THEME & LAYOUT (Optimized for Small Figure)
  theme_classic(base_size = 10) + 
  labs(
    x = "Parameter Estimate (95% HPD)",
    y = "Replicate (Sorted by Pruned Magnitude)",
    color = NULL 
  ) +
  
  theme(
    # Force White Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    # Clean Clutter
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", margin = margin(t = 3)),
    
    # Facet Headers
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10, hjust = 0),
    
    # Legend (Top)
    legend.position = "top", 
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    
    # Margins
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  )
