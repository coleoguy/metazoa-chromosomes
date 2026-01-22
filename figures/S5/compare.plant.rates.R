new_dir <- "../results/plants.rerun/"
old_dir <- "../results/exponential.prior - full model/mentor_results/"
ages    <- read.csv("../data/clade.ages.csv")
files   <- list.files(new_dir, pattern="\\.csv$", full.names=FALSE)

dys_log <- function(file, dir){
  clade <- strsplit(file, ".", fixed=TRUE)[[1]][1]
  age   <- ages$Age[match(tolower(clade), tolower(ages$Clade))]
  dat   <- read.csv(file.path(dir, file))
  dat   <- dat[dat$i > 100, , drop=FALSE]
  rate  <- rowSums(dat[, c("asc1","desc1")]) / age
  log10(rate[rate > 0 & is.finite(rate)])
}

par(mfrow=c(2,4), mar=c(4,4,3,1))

for(file in files){
  clade <- strsplit(file, ".", fixed=TRUE)[[1]][1]
  
  old_x <- dys_log(file, old_dir)
  new_x <- dys_log(file, new_dir)
  
  old_d <- density(old_x, cut=0.5)
  new_d <- density(new_x, cut=0.5)
  
  xlim <- range(old_d$x, new_d$x)
  ylim <- range(old_d$y, new_d$y)
  
  plot(old_d, xlim=xlim, ylim=ylim, main=clade, cex.main=1.5, cex.axis = 1.2,
       cex.lab = 1.2,
       xlab="Log dysploidy rate (events/Myr)", ylab="Density")
  lines(new_d, lty=2)
  abline(v=mean(old_x), lwd=2)
  abline(v=mean(new_x), lwd=2, lty=2)
  legend("topleft", legend=c("old","new"), lty=c(1,2), bty="n")
}
