### Megan Copeland
### Plotting chrom variance plants vs animals by dysploidy

# function to clean haploid column (to handle ranges, comma seperated values, multiple entries)
CleanHaploid <- function(df) {
  stoch_round <- function(x) floor(x) + (runif(1) < (x - floor(x)))
  for(i in 1:nrow(df)){
    x <- df$haploid[i]
    if(is.numeric(x)){
      if(!x == round(x)){
        df$haploid[i] <- stoch_round(x)
      }
    }else{
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

files <- list.files("../../data/chrome/", pattern = "\\.csv$", full.names = TRUE)
rates <- read.csv("../../results/median_dysploidy_rates.csv", stringsAsFactors = FALSE)

plants <- c("asteraceae", "brassicaceae", "fabaceae", "gymnospermae", "solanaceae",
            "liliaceae", "pteridophyta", "rubiaceae", "bryophyta", "orchidaceae",
            "passifloraceae")

# collect variance values
var_dat <- data.frame(
  Clade = character(),
  Variance = numeric(),
  Group = character())

for (i in seq_along(files)) {
  dat <- read.csv(files[i], stringsAsFactors = FALSE)
  dat <- CleanHaploid(dat)
  
  clade_name <- tools::file_path_sans_ext(basename(files[i]))
  clade_var  <- var(dat$haploid, na.rm = TRUE)
  clade_group <- if (clade_name %in% plants) "Plant" else "Animal"
  
  var_dat <- rbind(
    var_dat,
    data.frame(
      Clade = clade_name,
      Variance = clade_var,
      Group = clade_group))
}

# merge with dysploidy rates
rates$Clade <- as.character(rates$Clade)
plot_dat <- merge(var_dat, rates[, c("Clade", "Median_Rate")], by = "Clade")
plot_dat <- plot_dat[-36,]

plants_dat  <- plot_dat[plot_dat$Group == "Plant", ]
animals_dat <- plot_dat[plot_dat$Group == "Animal", ]

par(mfrow = c(1,2))   # 1 row, 2 plots

# Plants
plot(plants_dat$Median_Rate,
     plants_dat$Variance,
     pch = 16,
     xlab = "Median dysploidy rate",
     ylab = "Variance in haploid chromosome number",
     main = "Plants")

# Animals
plot(animals_dat$Median_Rate,
     animals_dat$Variance,
     pch = 16,
     xlab = "Median dysploidy rate",
     ylab = "Variance in haploid chromosome number",
     main = "Animals")


wilcox.test(Variance ~ Group, data = plot_dat)
