# Bandwidth Selection under Correlated Errors

This repository contains R code for the STAT\*6920 final project comparing bandwidth selection methods for local polynomial regression under correlated errors: Modified Cross-Validation (MCV), Partitioned Cross-Validation (PCV), and the Correlated Data Plug-In (CDPI) estimator.

## Repository Structure

```
├── src/
│   ├── functions.R       # Core implementations: bandwidth selectors, local polynomial fit, simulation driver
│   ├── sim.R             # Sequential simulation script (development/testing)
│   ├── sim_parallel.R    # Parallel simulation script (full B = 500 runs)
│   ├── visuals.R         # Result loading and figure generation
│   └── test.R            # Single-dataset sanity checks and Figure 1 (LOOCV breakdown plot)
├── data/                 # Saved simulation output (.rds files, not tracked)
├── plots                 # plots produced from `visuals.R` 
└── project.Rproj         # RStudio project file
```

## Core Implementations

All bandwidth selectors and fitting utilities are in `src/functions.R`:

- `r_true()` : True regression function r(x) = x + 2exp(−16x²)
- `generate_ar1_errors()` : AR(1) error generation with parameter ρ
- `local_poly_est()` : Local polynomial regression at a single point, returns coefficients and smoother row
- `bandwidth_cv()` : Ordinary LOOCV bandwidth selector
- `bandwidth_mcv()` : Modified cross-validation with block size parameter `l`
- `bandwidth_pcv()` : Partitioned cross-validation with partition parameter `g`, includes g^(−1/5) rescaling
- `bandwidth_cdpi()` : Correlated data plug-in selector using AR(1) pilot estimates and cubic derivative fit
- `simulate_bandwidth()` : Simulation driver: runs B replications in parallel across all methods and returns per-replication bandwidth and MISE matrices

## Running the Simulation

1. Run full parallel simulation: `source("src/sim_parallel.R")`  
   - Uses all available cores minus one, saves results to `data/` with datestamped filenames
2. Generate figures: `source("src/visuals.R")`  
   - Loads `.rds` files from `data/` and saves plots to `plots/`
3. Sequential version (slower, for debugging): `source("src/sim.R")`


## Simulation Design

- True regression function: r(x) = x + 2exp(−16x²)
- Error structure: AR(1) with ρ ∈ {0, 0.3, 0.6, 0.9}
- Sample sizes: n ∈ {100, 500, 1000}
- Variance: σ² = 0.5
- Replications: B = 500
- Bandwidth grid: H = {0.01, 0.011, ..., 0.20} (step 0.001)
- MCV block sizes: ℓ ∈ {1, 2, 5}
- PCV partition parameters: g ∈ {3, 5, 10}

## Dependencies

`MASS`, `parallel`, `pbapply`, `tidyverse`, `ggplot2`, `gridExtra`

## Authors

Vinay Joshy & Jackson Mitchell  
University of Guelph, Winter 2026