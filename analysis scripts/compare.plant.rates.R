### Megan Copeland
### Jan 20, 2026
### Comparing dysploidy rates on plant clades for trees before and after handling polytomies

# path to mcmc results from new run
new_dir <- "../results/plants.rerun/"

# path to mcmc results from old run
old_dir <- "../results/exponential.prior - full model/mentor_results/"

# path to clade ages 
ages    <- read.csv("../data/clade.ages.csv")

# get all files in the new run dir
files   <- list.files(new_dir, pattern="\\.csv$", full.names=FALSE)

# loop through all files 
for (i in 1:length(files)) {
  file <- files[[i]] # grab file i
  clade <- strsplit(file, ".", fixed=TRUE)[[1]][1] # strip everything but clade name
  age   <- ages$Age[match(tolower(clade), tolower(ages$Clade))] # get clade age
  
  old_dat <- read.csv(file.path(old_dir, file)) #read old mcmc file
  old_dat <- old_dat[old_dat$i > 100, , drop=FALSE] # discard burn in
  old_rate <- rowSums(old_dat[, c("asc1","desc1")]) / age # convert rates
  old_x <- log10(old_rate[old_rate > 0]) # log rates above zero
  
  new_dat <- read.csv(file.path(new_dir, file))
  new_dat <- new_dat[new_dat$i > 100, , drop=FALSE]
  new_rate <- rowSums(new_dat[, c("asc1","desc1")]) / age
  new_x <- log10(new_rate[new_rate > 0])
  
  old_d <- density(old_x, cut=0.5)
  new_d <- density(new_x, cut=0.5)
  
  plot(old_d, xlim=range(old_d$x, new_d$x), ylim=range(old_d$y, new_d$y),
       main=clade, xlab="log10(dysploidy rate; events/Myr)", ylab="Density")
  
  lines(new_d, lty=2)
  abline(v=mean(old_x), lwd=2)
  abline(v=mean(new_x), lwd=2, lty=2)
  legend("topright", legend=c("old","new"), lty=c(1,2), bty="n")
}
