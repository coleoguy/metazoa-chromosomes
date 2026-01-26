# ------------------------------------------------------------
# Full script: distinct color per clade, NO legend,
# and non-overlapping right-side colored clade labels with leader lines
# ------------------------------------------------------------

# ---- 0. Load data ----
animals <- read.csv("animals.csv")
plants  <- read.csv("plants.fungi.csv")
rates   <- read.csv("../../results/median_dysploidy_rates.csv")
rates$Clade <- paste0(toupper(substr(rates$Clade, 1, 1)), substr(rates$Clade, 2, nchar(rates$Clade)))

# ------------------------------------------------------------
# Full script: distinct color per clade, NO legend,
# and non-overlapping right-side colored clade labels with leader lines
# ------------------------------------------------------------

# ---- 1. Build merged plotting dataset ----
new_dat <- data.frame()

for (i in 1:length(rates$Clade)) {
  
  clade <- rates$Clade[i]
  rate  <- rates$Median_Rate[i]
  
  # ---- CASE 1: plants ----
  if (clade %in% plants$Family) {
    
    idx <- plants$Family == clade
    
    rows <- data.frame(
      genome_size = plants$Cvalue[idx],
      Clade       = clade,
      Median_Rate = rate,
      Source      = "Plants"
    )
    
    new_dat <- rbind(new_dat, rows)
    
  } else {
    
    # ---- CASE 2: animals ----
    family_name <- clade
    
    if (clade == "Iguania")    family_name <- "Iguanidae"
    if (clade == "Scincoidea") family_name <- "Scincidae"
    
    idx <- animals$Order  == family_name |
      animals$Family == family_name |
      animals$Class  == family_name
    
    if (any(idx, na.rm = TRUE)) {
      
      rows <- data.frame(
        genome_size = animals$Cvalue[idx],
        Clade       = clade,
        Median_Rate = rate,
        Source      = "Animals"
      )
      
      new_dat <- rbind(new_dat, rows)
    }
  }
}

# ---- 2. Libraries ----
library(ggplot2)
library(dplyr)

# ---- 3. Clean / coerce types ----
new_dat$genome_size <- as.numeric(new_dat$genome_size)
new_dat$Median_Rate <- as.numeric(new_dat$Median_Rate)

# Drop rows that break log scales
new_dat <- new_dat %>%
  filter(
    is.finite(genome_size), is.finite(Median_Rate),
    genome_size > 0, Median_Rate > 0
  )

# Stable factor ordering
new_dat$Clade <- factor(new_dat$Clade, levels = sort(unique(new_dat$Clade)))
n_clades <- nlevels(new_dat$Clade)

# ---- 4. One unique color per clade ----
cols <- setNames(
  grDevices::hcl.colors(n_clades, "Dark 3"),
  levels(new_dat$Clade)
)

# ---- 5. Build non-overlapping label positions + leader lines ----
# Representative y per clade: median rate
# Representative x per clade: a high quantile of genome_size (anchor near right edge of that clade’s cloud)
rep_df <- new_dat %>%
  group_by(Clade) %>%
  summarize(
    y_anchor = median(Median_Rate, na.rm = TRUE),
    x_anchor = as.numeric(quantile(genome_size, probs = 0.90, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  filter(is.finite(x_anchor), is.finite(y_anchor), x_anchor > 0, y_anchor > 0)

# Global plot bounds
x_max <- max(new_dat$genome_size, na.rm = TRUE)
y_min <- min(new_dat$Median_Rate, na.rm = TRUE)
y_max <- max(new_dat$Median_Rate, na.rm = TRUE)

# Put labels in a neat vertical stack *in log space* (so it looks even on a log axis)
# Order by the clade’s anchor y so leader lines cross less.
rep_df <- rep_df %>%
  arrange(y_anchor) %>%
  mutate(
    # label column x location slightly beyond the data
    x_label = x_max * 1.25,
    # evenly spaced positions in log10(y)
    y_label = 10^(seq(log10(y_min), log10(y_max), length.out = n()))
  )

# ---- 6. Plot ----
p <- ggplot(new_dat, aes(x = genome_size, y = Median_Rate, color = Clade)) +
  geom_point(alpha = 0.3, size = 1.8, show.legend = FALSE) +
  
  # leader lines from clade anchor -> label position
  geom_segment(
    data = rep_df,
    aes(x = x_anchor, y = y_anchor, xend = x_label, yend = y_label, color = Clade),
    linewidth = 0.4,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  # stacked colored labels
  geom_text(
    data = rep_df,
    aes(x = x_label, y = y_label, label = Clade, color = Clade),
    hjust = 0,
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_x_log10(expand = expansion(mult = c(0.02, 0.35))) +
  scale_y_log10() +
  scale_color_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "C-Value (pg)",
    y = "Dysploidy rate (Median_Rate)"
  )

print(p)
