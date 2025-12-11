library(ape)
library(diversitree)
library(chromePlus)
library(phangorn)


# Simulate a random tree (20 tips)
# use compute.brlen to ensure it is ultrametric (required for diversitree)
n_tips <- 20
tree <- rtree(n = n_tips)
tree <- compute.brlen(tree, method = "Grafen") 

# Simulate random chromosome data
# create a dataframe matching tips to random haploid numbers (e.g., 10 to 20)
dat <- data.frame(
  species = tree$tip.label,
  haploid = sample(10:20, n_tips, replace = TRUE))


# Convert the chromosome data to ChromePlus format matrix
# 'buffer' allows for chromosome numbers slightly outside range
rng <- c(min(dat$haploid) - 1, max(dat$haploid) + 1)
mat <- datatoMatrix(x = dat, buffer = 1, hyper = FALSE)

# overwrite the actual data matrix with 1s.
mat[,] <- 1/ncol(mat)
mat <- as.matrix(mat)


lik <- make.mkn(tree = tree, states = mat, k = ncol(mat), 
                strict = FALSE,
                control = list(method = "ode", root = ROOT.OBS))

# Constrain the likelihood
conlik <- constrainMkn(data = mat, lik = lik, hyper = FALSE,
                       polyploidy = FALSE, verbose = FALSE)


# Running a short chain for demonstration
iter <- 1000 

res <- mcmc(lik = conlik,
            x.init = runif(length(argnames(conlik))), # random start values
            prior = make.prior.exponential(2),        # exponential prior
            nsteps = iter, 
            w = 1)

## ----- RESULTS -----
# View the first few lines of the result
print(head(res))

# Plot the trace of the first parameter to check if it moves
plot(res$asc1, type='l', main="Trace of asc1 (Prior Only)", ylab="Rate", xlab="Generation")

# If you want to save it:
# write.csv(res, "simulated_prior_run.csv", row.names = FALSE)