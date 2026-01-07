library(coda)

## ---- user inputs ----
dir1 <- "../results/exponential.prior - full model/mentor_results/"
dir2 <- "../results/uniform prior - full model/"
ages <- read.csv("../data/clade.ages.csv")

## ---- get file lists ----
file.list1 <- sort(list.files(dir1, full.names=TRUE))
file.list2 <- sort(list.files(dir2, full.names=TRUE))

## ---- results data frame ----
out <- data.frame(
  clade = ages$Clade,
  exp_hpd_low = NA_real_,
  exp_hpd_high = NA_real_,
  exp_median = NA_real_,
  uni_hpd_low = NA_real_,
  uni_hpd_high = NA_real_,
  uni_median = NA_real_,
  stringsAsFactors = FALSE
)

## ---- loop in requested form ----
for (i in 1:length(file.list1)) {
  
  exp.res <- rowSums(read.csv(file.list1[i], stringsAsFactors = FALSE, check.names = FALSE)[,2:3])/ages$Age[i]
  h1 <- HPDinterval(as.mcmc(exp.res))
  out$exp_hpd_low[i]  <- h1[1, "lower"]
  out$exp_hpd_high[i] <- h1[1, "upper"]
  out$exp_median[i]   <- median(exp.res)
  
  uni.res <- rowSums(read.csv(file.list2[i], stringsAsFactors = FALSE, check.names = FALSE)[,2:3])/ages$Age[i]
  h2 <- HPDinterval(as.mcmc(uni.res))
  out$uni_hpd_low[i]  <- h2[1, "lower"]
  out$uni_hpd_high[i] <- h2[1, "upper"]
  out$uni_median[i]   <- median(uni.res)
}



res <- out

x <- res[[4]]  # exp median
y <- res[[7]]  # unif median

xticks <- 10^pretty(log10(range(x)))
yticks <- 10^pretty(log10(range(y)))

fmt <- function(z) format(signif(z, 3), trim = TRUE, scientific = TRUE)

op <- par(no.readonly = TRUE)
on.exit(par(op))

## more left margin, and move y-axis title further left
par(mar = c(5, 8, 4, 2) + 0.1)      # bottom, left, top, right
par(mgp = c(2.6, 1.0, 0))           # title, labels, line (distances)
par(cex.axis = 0.9)

plot(
  x, y,
  log = "xy",
  xaxt = "n", yaxt = "n",
  pch = 16,
  xlab = "Median rate (events per My) — exponential prior",
  ylab = ""                         # add with mtext so we can control spacing
)

axis(1, at = xticks, labels = fmt(xticks))
axis(2, at = yticks, labels = fmt(yticks), las = 1)

mtext("Median rate (events per My) — uniform prior", side = 2, line = 5.5)

abline(a = 0, b = 1, lty = 2)