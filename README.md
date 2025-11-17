---
editor_options: 
  markdown: 
    wrap: 72
---

# CURE-chromosome-number

The in this repository recapitulate analyss of CURE chromosome Evolution
Project. This project investigates rate variation in fusion, fission,
and whole genome duplication.

1.  README.md: The file that you are reading
2.  **data**
3.  chrome: contains csv files with phenotype data for all clades
4.  trees: contains phylogenies for all clades
5.  **figures**: contains scripts to produce each figure.
6.  **script**: contains all scripts for data validation and analysis

TODO:

student scripts per clade

faculty script that runs all clades

supplemental file that has paragraph for each clade and plots of data

see about splitting pteridophytes

figure out whether it is correct to compare rate estimates across two
trees one that is small and one that is large that have the same prior.
Essentially could do this with a small sim:

simulate 100 phylogenies with chromsome data with 50 tips and a unit
length, use a symmetric transition matrix with a rate of 0.2, do this a
second time with a 200 tip tree. Use the same prior on both and ask if
we are biased or if it is just more noise on the small tree. My concern
is that a small tree the new data will not have the same amount of power
to pull the posterior away from the prior.
