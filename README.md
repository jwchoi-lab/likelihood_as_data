# Robust Model Selection using Likelihood as Data

This repository contains all code and data needed to reproduce the analyses, figures, and tables in the manuscript *Robust Model Selection using Likelihood as Data*. 

---


## Software requirements

- **R (≥ 4.3)** is required. We use `renv` to capture the exact package environment. To install all required packages, run from the repository root:

```r
Rscript -e "install.packages('renv'); renv::restore()"
```

- **Julia (≥ 1.9)** — optional; needed only to re-run the MFM posterior sampler (see below).
- **STRUCTURE 2.3.4** — optional; needed only to re-run the admixture analysis (see below).

---

## Reproducing the analyses

All scripts source `code/functions.R`, which contains the core LaD functions. Run each script from the repository root in R or RStudio.

| Script | Paper section | Output |
|--------|---------------|--------|
| `code/examples/01_shapley_gmm.R`          | Section 5.1      | Figures 1, 3 |
| `code/examples/02_sparseMVN.R`            | Section 5.2      | Figures 4, 5 |
| `code/examples/03_tpc.R`                  | Section 5.3      | Figure 6     |
| `code/examples/04_structure_admixture.R`  | Section S3       | Figure S1, S2    |
| `code/theory/instability.R`               | Section 4.2      | Figure 2     |

All figures are saved to `output/figures`.


### Reproducibility notes

Precomputed results are provided in `output/` for scripts that involve random initialization or Monte Carlo sampling, so loading these files produces figures identical to those in the paper. Results may differ slightly if re-run from scratch, but qualitative conclusions are unchanged. 

- **`01_shapley_gmm.R`**: Precomputed fits in `output/shapley_gmm_fitted.Rdata` (Gaussian mixture EM across sample sizes). The MFM posterior uses precomputed Julia output from `data/processed/julia_run/`.
- **`02_sparseMVN.R`**: Precomputed results in `output/sparseMVN_fitted.Rdata`.
- **`03_tpc.R`**: Precomputed model fits in `output/tpc_fitted.Rdata`.
- **`04_structure_admixture.R`**: Uses precomputed results in `output/admixture_fitted.Rdata` and uses STRUCTURE output in `data/processed/structure_run/`.

---

### Julia (optional)

The MFM posterior was computed using [BayesianMixtures.jl](https://github.com/jwmi/BayesianMixtures.jl) (Miller & Harrison 2018). Precomputed summaries are in `data/processed/julia_run/`, so Julia is not required to reproduce the figures. To re-run the sampler:

1. Install Julia from [julialang.org](https://julialang.org/downloads/)
2. Install BayesianMixtures.jl:
   ```julia
   using Pkg; Pkg.add(url="https://github.com/jwmi/BayesianMixtures.jl")
   ```
3. See `code/julia/` for the sampler script and environment files.

---


### STRUCTURE (optional)

The population structure admixture models (Section S3) uses STRUCTURE 2.3.4 to fit the admixture model. All STRUCTURE outputs needed for the figures are already included in `data/processed/structure_run/`, so STRUCTURE is not required to reproduce the results. 
The instructions below are only for readers who wish to re-fit the model.

**Software**
- STRUCTURE 2.3.4 — download from [pritchardlab.stanford.edu](https://web.stanford.edu/group/pritchardlab/structure.html)
- Binary location: `./structure` (or on PATH)

**Inputs**
- Parameter files: `code/structure/mainparams`, `code/structure/extraparams`
- Genotype file: `data/raw/admixture/brooktrout.txt` (ensure `mainparams` points to this file)

**Outputs**
- Written to: `data/processed/structure_run/`
- Naming convention: `results_K{K}_rep{rep}` (plus STRUCTURE suffixes)

**Example command-line loop**
```bash
mkdir -p data/processed/structure_run

for K in {1..10}; do
  for rep in {1..20}; do
    SEED=$((1000 * K + rep))
    ./structure \
      -K $K \
      -D $SEED \
      -m code/structure/mainparams \
      -e code/structure/extraparams \
      -o data/processed/structure_run/results_K${K}_rep${rep}
  done
done
```

---


## Data sources

All datasets are included to facilitate reproduction of the results. Please cite the original sources if you reuse these data, and respect their licenses.

| Dataset | File | Source | License |
|---------|------|--------|---------|
| Shapley galaxy radial velocities | `data/raw/shapley/Shapley_galaxy.dat` | Drinkwater et al. (2004), data: [Penn State Astrostatistics](https://sites.psu.edu/astrostatistics/datasets-shapley-galaxy-dataset/) | No restrictions; please cite source |
| *P. minimum* growth rates | `data/raw/tpc/thermal_performance_datasets.csv` | Kontopoulos et al. (2024), Figshare [doi:10.6084/m9.figshare.24106161.v3](https://doi.org/10.6084/m9.figshare.24106161.v3) | **CC BY 4.0** |
| Brook trout microsatellite genotypes | `data/raw/admixture/brooktrout.txt` | Erdman et al. (2022); data: [USGS ScienceBase](https://www.sciencebase.gov/catalog/item/611d264cd34e40dd9c01284e) | **U.S. public domain** (USGS data release) |

---



## References for software and data

- R Core Team (2021). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing.
- Scrucca et al. (2016). mclust 5. *The R Journal*, 8(1), 289.
- Padfield et al. (2021). rTPC and nls.multstart. *Methods in Ecology and Evolution*, 12(6), 1138–1143.
- Miller & Harrison (2018). Mixture models with a prior on the number of components. *Journal of American Statistical Association*, 113(521), 340–356.
- Pritchard et al. (2000). Inference of population structure using multilocus genotype data. *Genetics*, 155(2), 945–959.
- Drinkwater et al. (2004). The large scale distribution of galaxies in the Shapley supercluster. *Publications of the Astronomical Society of Australia*, 21(1), 89–96.
- Kontopoulos et al. (2024). No universal mathematical model for thermal performance curves. *Nature Communications*, 15(1), 8855. 
- Grzebyk & Berland (1996). Influences of temperature, salinity and irradiance on growth of *Prorocentrum minimum*. *Journal of Plankton Research*, 18(10), 1837–1849.
- Erdman et al. (2022). Broadscale population structure and hatchery introgression of Midwestern Brook Trout. *Transactions of the American Fisheries Society*, 151(1), 81–99.
