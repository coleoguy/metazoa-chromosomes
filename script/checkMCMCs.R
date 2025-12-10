library(coda)
files <- list.files("../results",, full.names=T)
for(i in 1:length(files)){
  dat <- read.csv(files[i])
  effectiveSize(as.mcmc(dat))
}