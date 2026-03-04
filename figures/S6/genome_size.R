library(ggplot2)

# ---- Load data ----
animals <- read.csv("animals.csv")
plants  <- read.csv("plants.fungi.csv")
rates   <- read.csv("../../results/median_dysploidy_rates.csv")

rates$Clade <- paste0(
  toupper(substr(rates$Clade, 1, 1)),
  substr(rates$Clade, 2, nchar(rates$Clade))
)

# ---- Build dataset ----
new_dat <- data.frame(stringsAsFactors = FALSE)

for (i in seq_len(nrow(rates))) {
  
  clade <- rates$Clade[i]
  rate  <- rates$Median_Rate[i]
  
  if (clade %in% plants$Family) {
    
    idx <- plants$Family == clade
    new_dat <- rbind(new_dat, data.frame(
      genome_size = plants$Cvalue[idx],
      Clade       = clade,
      Median_Rate = rate,
      stringsAsFactors = FALSE
    ))
    
  } else {
    
    family_name <- clade
    
    idx <- animals$Order  == family_name |
      animals$Family == family_name |
      animals$Class  == family_name
    
    if (any(idx, na.rm = TRUE)) {
      new_dat <- rbind(new_dat, data.frame(
        genome_size = animals$Cvalue[idx],
        Clade       = clade,
        Median_Rate = rate,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# ---- Clean ----
new_dat$genome_size <- as.numeric(new_dat$genome_size)
new_dat$Median_Rate <- as.numeric(new_dat$Median_Rate)

# Convert pg -> Mb (1 pg ≈ 978 Mbp)
new_dat$genome_size <- new_dat$genome_size * 978

new_dat$Clade <- factor(new_dat$Clade, levels = sort(unique(new_dat$Clade)))

# ---- Anchor: max genome_size per clade; y = clade's (constant) Median_Rate ----
rep_df <- aggregate(genome_size ~ Clade, data = new_dat, FUN = max)
names(rep_df)[names(rep_df) == "genome_size"] <- "x_anchor"
rep_df$y_anchor <- new_dat$Median_Rate[match(rep_df$Clade, new_dat$Clade)]

# ---- Label positions (stacked evenly in log space) ----
x_max <- max(new_dat$genome_size, na.rm = TRUE)
y_min <- min(new_dat$Median_Rate, na.rm = TRUE)
y_max <- max(new_dat$Median_Rate, na.rm = TRUE)

rep_df <- rep_df[order(rep_df$y_anchor), , drop = FALSE]
rep_df$x_label <- x_max * 1.30
rep_df$y_label <- 10^(seq(log10(y_min), log10(y_max), length.out = nrow(rep_df)))

# ---- Colors (fixed per clade) ----
n_clades <- nlevels(new_dat$Clade)
cols <- setNames(grDevices::hcl.colors(n_clades, "Dark 3"), levels(new_dat$Clade))

# ---- Plot ----
p <- ggplot(new_dat, aes(x = genome_size, y = Median_Rate, color = Clade)) +
  geom_point(alpha = 0.35, size = 1.8, show.legend = FALSE) +
  
  geom_segment(
    data = rep_df,
    aes(x = x_anchor, y = y_anchor, xend = x_label, yend = y_label, color = Clade),
    inherit.aes = FALSE,
    linewidth = 0.4,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  
  geom_text(
    data = rep_df,
    aes(x = x_label, y = y_label, label = Clade, color = Clade),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3,
    show.legend = FALSE
  ) +
  
  scale_x_log10(expand = expansion(mult = c(0.02, 0.35))) +
  scale_y_log10() +
  scale_color_manual(values = cols) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Genome size (Mb)", y = "Dysploidy rate (Median_Rate)")

print(p)

### ---- PGLS ----
library(ape)
library(nlme)

tree <- read.tree("cladetree.new")

# One row per clade: mean genome size + clade dysploidy rate
dat_clade <- data.frame(
  Clade = unique(new_dat$Clade),
  Median_Rate = new_dat$Median_Rate[match(unique(new_dat$Clade), new_dat$Clade)],
  GenomeSize = as.numeric(tapply(new_dat$genome_size, new_dat$Clade, mean, na.rm = TRUE)),
  stringsAsFactors = FALSE)

# Match tree and data
rownames(dat_clade) <- dat_clade$Clade
shared <- intersect(tree$tip.label, dat_clade$Clade)

tree_p <- drop.tip(tree, setdiff(tree$tip.label, shared))
dat_p  <- dat_clade[tree_p$tip.label, , drop = FALSE]

# Fit PGLS (Brownian motion correlation)
phylo_cor <- corBrownian(1, phy = tree_p)

pgls_model <- gls(
  Median_Rate ~ GenomeSize,
  data = dat_p,
  correlation = phylo_cor,
  method = "ML")

summary(pgls_model)


dat_p$pred <- predict(pgls_model)
dat_plot <- dat_p[order(dat_p$GenomeSize), ]

ggplot(dat_plot, aes(GenomeSize, Median_Rate)) +
  geom_point(size = 2) +
  geom_line(aes(y = pred), color = "red", linewidth = 1) +
  scale_y_log10() +          # <-- match your original plot
  theme_classic()
