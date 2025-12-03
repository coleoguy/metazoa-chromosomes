# ============================
# Phylogeny + Beeswarm Points + Clade Labels (colored by higher classification)
# ============================

# ---- 1. Libraries ----
library(ape)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(dplyr)
library(purrr)
library(stringr)
library(viridisLite)
library(readr)

# ---- 2. Load haploid chromosome data ----
file_list <- list.files(path = "../../data/chrome", pattern = "\\.csv$", full.names = TRUE)

all_data <- list()
for (file in file_list) {
  df <- read.csv(file)
  if (!"haploid" %in% names(df)) {
    cat("Missing 'haploid' column in:", basename(file), "\n")
    next
  }
  file_name <- tools::file_path_sans_ext(basename(file))
  temp_df <- data.frame(clade = tolower(file_name), haploid = df$haploid)
  all_data[[length(all_data) + 1]] <- temp_df
}
dat <- do.call(rbind, all_data)

# ---- 3. Load tree ----
tree <- read.tree("cladetree.new")
tree$tip.label <- tolower(tree$tip.label)

# ---- 4. Clean haploid values ----
dat_clean <- dat %>%
  mutate(
    haploid = str_trim(haploid),
    haploid = map_dbl(haploid, ~ {
      vals <- str_split(.x, ",|–|-")[[1]] |> str_trim()
      nums <- suppressWarnings(as.numeric(vals))
      nums <- nums[!is.na(nums) & nums > 0]
      if (length(nums) == 0) return(NA_real_)
      if (length(nums) == 1) return(nums)
      if (length(nums) == 2 && diff(nums) > 0) return(runif(1, nums[1], nums[2]))
      return(sample(nums, 1))
    })
  ) %>%
  filter(is.finite(haploid), haploid > 0)

# ---- 5. Match data to tree tips ----
final_df <- dat_clean %>%
  filter(clade %in% tree$tip.label) %>%
  mutate(
    clade = factor(clade, levels = tree$tip.label),
    log_haploid = log10(haploid)
  )

# ---- 5b. Add higher classification from CSV ----
class_df <- read.csv("higher_class.csv", stringsAsFactors = FALSE)
class_df <- class_df %>% mutate(Clade = tolower(str_trim(Clade)))

final_df <- final_df %>%
  mutate(clade = tolower(str_trim(clade))) %>%
  left_join(class_df, by = c("clade" = "Clade")) %>%
  rename(higher_class = Higher.Classification)

# ---- 6. Calculate species counts ----
counts <- final_df %>%
  count(clade, name = "count") %>%
  mutate(label_text = paste0(clade, " (n=", count, ")"))

# ---- 7. Build base tree ----
p <- ggtree(tree) +
  theme_tree() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# ---- 8. Align y positions with actual plotted order ----
tip_order <- p$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

baseline_df <- data.frame(clade = tree$tip.label) %>%
  mutate(y_pos = match(clade, tip_order))

final_df <- final_df %>%
  mutate(y_pos = match(as.character(clade), tip_order))

counts <- counts %>%
  mutate(y_pos = match(as.character(clade), tip_order))

# ---- 9. Get x range ----
x_min <- min(final_df$log_haploid)
x_max <- max(final_df$log_haploid) + 0.25

# ---- 10. Final plot ----
p_final <- p +
  geom_facet(
    data = baseline_df,
    panel = "Haploid counts",
    mapping = aes(y = y_pos, yend = y_pos, x = x_min, xend = x_max),
    geom = geom_segment,
    linewidth = 0.35,
    color = "grey55",
    alpha = 0.6
  ) +
  geom_facet(
    data = final_df,
    panel = "Haploid counts",
    mapping = aes(x = log_haploid, y = y_pos, color = higher_class),
    geom = geom_point,
    position = position_jitter(width = 0.2, height = 0.2),
    size = 1.2,
    alpha = 0.1
  ) +
  geom_facet(
    data = counts %>% mutate(x_text = 0),
    panel = "Clade labels",
    mapping = aes(x = x_text, y = y_pos, label = label_text),
    geom = geom_text,
    size = 2.6,
    hjust = 0
  ) +
  scale_color_manual(
    name = "Higher classification",
    values = c(
      "Chondrichthyes" = "#414487",
      "Actinopterygii" = "#2A788E",
      "Amphibia" = "#22A884",
      "Mammalia" = "#7AD151",
      "Reptilia" = "#F8961E",
      "Arachnida" = "#D1495B",
      "Insecta" = "#8E3B9E",
      "Angiosperm" = "#1F968B",
      "Gymnosperms" = "#89C2D9",
      "Pteridophytes" = "#5E4FA2",
      "Bryophyta" = "#277DA1"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 3)
    )
  ) +
  scale_x_continuous(
    breaks = log10(c(2, 5, 10, 20, 50, 100, 200)),
    labels = c(2, 5, 10, 20, 50, 100, 200)
  ) +
  labs(x = "Haploid chromosome number (log10)") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 12, margin = margin(t = 6)),
    axis.text.x = element_text(size = 10)
  ) +
  expand_limits(x = x_max + 0.2)

# ---- 11. Render ----
print(p_final)
