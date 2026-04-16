# Dismantling Chromosomal Stasis Across the Eukaryotic Tree of Life

Data, phylogenies, and analysis code for Copeland et al. (2026), which estimates rate variation in fusion, fission, polyploidy, and demi-polyploidy across 56 clades spanning the eukaryotic tree of life.

Karyotype records were collected by undergraduate researchers through Course-based Undergraduate Research Experiences (CUREs) at Texas A&M University. An interactive browser for the full dataset is available at [coleoguy.github.io/cures-karyotype-database.html](https://coleoguy.github.io/cures-karyotype-database.html).

## Repository layout

```
metazoa-chromosomes/
├── data/
│   ├── chrome/             # per-clade CSVs of species × haploid number
│   ├── trees/              # per-clade source phylogenies (Newick/Nexus)
│   ├── pruned-trees/       # trees pruned to match chromosome-data tips
│   ├── clade.ages.csv      # crown ages used in rate calculations
│   ├── tracking-sheet.xlsx # CURES student dataset provenance
│   └── cures-karyotype-database.csv   # full 63,542-row combined dataset
├── analysis-scripts/       # core MCMC rate estimation and diagnostics
├── adequacy-scripts/       # sampling and tree-size sensitivity analyses
├── results/                # MCMC output, AIC tables, posterior predictive sims
└── figures/                # scripts and source files for each manuscript figure
```

## Data

`data/cures-karyotype-database.csv` holds the complete assembled dataset: 63,542 species × haploid number records across 56 clades, with one karyotype-source citation per clade.

| column | description |
|---|---|
| `clade` | clade name (e.g. *Accipitriformes*, *Fabaceae*) |
| `species` | species binomial |
| `haploid_number` | haploid chromosome number (n) |
| `citation` | karyotype source for the clade |

Per-clade CSVs in `data/chrome/` are the analysis-ready files used by the R scripts.

## Reproducing the analyses

Primary dependencies: R (≥ 4.2), `chromePlus`, `diversitree`, `ape`, `phytools`, `coda`.

From `analysis-scripts/`:

1. `produce_pruned_trees.r` — match chromosome datasets to trees, write pruned trees.
2. `parallele.analysis.exponential.r` and `parallele.analysis.uniform.r` — run the full chromosome-evolution model across all clades under each prior.
3. `ESS.checker.R`, `checkMCMCs.R` — MCMC diagnostics.
4. `generate.PPS.R` — posterior predictive simulations per clade.
5. `model.choice.sensitivity.R`, `best.model.R`, `uni-vs-exp.R` — model comparison across parameterizations and priors.

Adequacy checks (`adequacy-scripts/`) assess how tree size and sampling influence posterior rate estimates using the scarab dataset and matched simulations.

## Citation

If you use this dataset or any of these scripts, please cite:

Copeland, M., McConnell, M., Barboza, A., Abraham, H.M., Alfieri, J., Arackal, S., Bernard, C.E., Bryant, K., Cast, S., Chien, S., Clark, E., Cruz, C.E., Diaz, A.Y., Deiterman, O., Girish, R., Harper, K., Hjelmen, C.E., Thompson, M.J., Koehl, R., Koneru, T., Laird, K., Lee, Y., Lopez, V.R., Murphy, M., Perez, N., Schmalz, S., Sylvester, T., and Blackmon, H. (2026). Dismantling Chromosomal Stasis Across the Eukaryotic Tree of Life. *bioRxiv* 2026.04.14.718287. https://doi.org/10.64898/2026.04.14.718287

```bibtex
@article{Copeland2026.04.14.718287,
  author    = {Copeland, Megan and McConnell, Meghann and Barboza, Andres and Abraham, Hannah M and Alfieri, James and Arackal, Steven and Bernard, Carrie E and Bryant, Kiedon and Cast, Shelbie and Chien, Sean and Clark, Emily and Cruz, Cassandra E and Diaz, Aileen Y and Deiterman, Olivia and Girish, Riya and Harper, Kaya and Hjelmen, Carl E and Thompson, Michelle J and Koehl, Rachel and Koneru, Tanvi and Laird, Kenzie and Lee, Yoonseo and Lopez, Virginia R and Murphy, Mallory and Perez, Nayeli and Schmalz, Sarah and Sylvester, Terrence and Blackmon, Heath},
  title     = {Dismantling Chromosomal Stasis Across the Eukaryotic Tree of Life},
  year      = {2026},
  doi       = {10.64898/2026.04.14.718287},
  publisher = {Cold Spring Harbor Laboratory},
  journal   = {bioRxiv},
  elocation-id = {2026.04.14.718287},
  url       = {https://www.biorxiv.org/content/early/2026/04/16/2026.04.14.718287}
}
```

Please also cite the original karyotype sources listed in the `citation` column of `data/cures-karyotype-database.csv`.

## License

Code is released under the MIT License (see `LICENSE`). The assembled dataset is released under CC-BY-4.0; please cite the paper above when reusing the data.
