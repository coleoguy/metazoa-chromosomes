These scripts perform the core analyses and checks for estimating chromosome evolution rates across clades using MCMC-based Mk models.

final-fitmodel.R: Runs the primary chromosome evolution model for a single clade or across multiple trees, including data cleaning, tree pruning, ultrametric adjustment, MCMC estimation, and writing posterior samples to file.

parallele.analysis.exponential.r: Executes the full chromosome evolution model across all clades in parallel using an exponential prior on rate parameters, writing per-clade MCMC results to disk.

parallele.analysis.uniform.r: Executes the full chromosome evolution model across all clades in parallel using a uniform prior on rate parameters, allowing sensitivity comparisons with the exponential prior.

prior_only.R: Runs a prior-only MCMC on simulated data to confirm that parameters move and to visualize the behavior of the prior independent of empirical data.

ESS.checker.R: Calculates effective sample sizes (ESS) for key rate parameters across MCMC output files to assess chain mixing and convergence.

checkMCMCs.R: Performs additional MCMC diagnostics and summaries on posterior samples, including ESS calculations using the coda package.

generate.PPS.R: Conducts posterior predictive simulations for each clade by simulating chromosome data from posterior samples and comparing empirical variance and entropy to simulated expectations.

model.choice.sensitivity.R: Compares alternative model parameterizations by summarizing posterior rates, calculating medians, and ranking clades to assess sensitivity to model structure.

clade.paragraphs.R: Compiles clade-level descriptive statistics (taxonomy, sampling overlap, tree resolution, and ages) to support manuscript text and results summaries.

cladetree.new: Newick-formatted clade-level tree used for ordering analyses and generating summary statistics across clades.

produce_pruned_trees.R: Generates pruned phylogenies for each clade by matching chromosome datasets to their corresponding trees and writing pruned trees to disk for downstream analyses.

S-pull of the prior - emp.R: Uses empirical data to evaluate the impact of tree size on posterior rate estimates by fitting the full dataset and repeatedly refitting models after randomly pruning half of the tips.

S-pull of the prior.R: Uses simulated trees and chromosome data to assess how prior choice and reduced sampling affect posterior rate estimates by comparing full and pruned trees across repeated simulations.
