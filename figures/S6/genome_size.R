animals <- read.csv("animals.csv")
plants  <- read.csv("plants.fungi.csv")
rates   <- read.csv("rates.csv")

new_dat <- data.frame()  # will become: genome_size, Clade, Median_Rate, Source

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
  }
  
  # ---- CASE 2: animals ----
  else {
    
    family_name <- clade
    if (clade == "Iguania") family_name <- "Iguanidae"
    if (clade == "Scincoidea") family_name <- "Scincidae"
    
    idx <- animals$Order  == clade |
      animals$Family == family_name |
      animals$Class  == clade
    
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

library(ggplot2)

new_dat$genome_size <- as.numeric(new_dat$genome_size)


# make sure clade is a factor so legend only includes what’s present
new_dat$Clade <- factor(new_dat$Clade, levels = sort(unique(new_dat$Clade)))

# Plot (all points)
p <- ggplot(new_dat, aes(x = genome_size, y = Median_Rate, color = Clade)) +
  geom_point(alpha = 0.55, size = 1.8) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), ncol = 1)) +
  theme(
    legend.key.height = grid::unit(0.25, "lines"),
    legend.spacing.y  = grid::unit(0.05, "lines")
  ) +
  labs(
    x = "Genome size (Cvalue)",
    y = "Dysploidy rate (Median_Rate)",
    color = "Clade"
  )

print(p)
