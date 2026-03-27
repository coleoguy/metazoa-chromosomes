new_dir <- "../results/plants.rerun/"
old_dir <- "../results/exponential.prior - full model/mentor_results/"
ages    <- read.csv("../data/clade.ages.csv")
files   <- list.files(new_dir, pattern="\\.csv$", full.names=FALSE)

dys_rate <- function(file, dir){
  clade <- strsplit(file, ".", fixed=TRUE)[[1]][1]
  age   <- ages$Age[match(tolower(clade), tolower(ages$Clade))]
  dat   <- read.csv(file.path(dir, file))
  dat   <- dat[dat$i > 100, , drop=FALSE]
  rate  <- rowSums(dat[, c("asc1","desc1")]) / age
  rate[rate > 0 & is.finite(rate)]
}

# 8 plots + 1 legend panel (bottom right)
layout(matrix(c(1,2,3,4,
                5,6,7,8,
                9,9,9,9), nrow=3, byrow=TRUE),
       heights=c(1,1,0.25))  # smaller legend row

par(mar=c(4,4,3,1))

for(i in seq_along(files)){
  file  <- files[i]
  clade <- strsplit(file, ".", fixed=TRUE)[[1]][1]
  
  old_x <- dys_rate(file, old_dir)
  new_x <- dys_rate(file, new_dir)
  
  old_d <- density(old_x, cut=0.5)
  new_d <- density(new_x, cut=0.5)
  
  xlim <- range(old_d$x, new_d$x)
  ylim <- range(old_d$y, new_d$y)
  
  plot(old_d, xlim=xlim, ylim=ylim, main=clade,
       cex.main=1.5, cex.axis=1.2, cex.lab=1.2,
       xlab="Dysploidy rate (events/Myr)", ylab="Density")
  lines(new_d, lty=2)
  abline(v=mean(old_x), lwd=2)
  abline(v=mean(new_x), lwd=2, lty=2)
}

# ---- dedicated legend panel ----
par(mar=c(0,0,0,0))
plot.new()

legend("center",
       legend=c("unpruned","pruned"),
       lty=c(1,2),
       lwd=2,
       horiz=TRUE,
       bty="n",
       cex=1.2)