# Review of post-clustering inference methods

This repository contains scripts to reproduce all numerical experiments and figures from the paper:

> Enjalbert-Courrech, N., Maugis-Rabusseau, C., & Neuvial, P. (2022). *Review of Post-Clustering Inference Methods*. 

A PDF version of each figure is stored in the `figures/` directory, and the corresponding R scripts are available in `scripts/`.

## Generate all figures and numerical results

To generate all figures and numerical results at once, run:

```sh
chmod +x run_all_scripts.sh
./run_all_scripts.sh
```

## Run a single script

To run an individual script, use: 

```sh
R --vanilla < [script_name].R
```

## Install dependencies

To install the required R packages, ensure you have [`pak`](https://pak.r-lib.org/) installed, then run:

```r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# Load package names from requirements.txt
packages <- readLines("requirements.txt")

# Install all required packages
pak::pkg_install(packages)
```

