
library(coda)
coda_samples <- mcmc(samples[, -1], start = 1) 

# 3. Calculate Effective Sample Size
ess_values <- effectiveSize(coda_samples)

# View the results
print(ess_values)

# 3. Calculate Effective Sample Size

# View the results
print(ess_values)

