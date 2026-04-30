# bHIVE

## B-cell Hybrid Immune Variant Engine

<img align="right" src="man/figures/logo.png" width="305" height="352">

<!-- badges: start -->
[![R-CMD-check](https://github.com/BorchLab/bHIVE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BorchLab/bHIVE/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/BorchLab/bHIVE/graph/badge.svg)](https://app.codecov.io/gh/BorchLab/bHIVE)
<!-- badges: end -->

### Overview

**bHIVE** is an R package implementing a modular Artificial Immune System (AIS) framework for clustering and classification. Built on [AI-Net](https://www.dca.fee.unicamp.br/~vonzuben/research/lnunes_dout/artigos/DMHA.PDF) (de Castro & Von Zuben 2001), bHIVE extends the classical algorithm with biologically-grounded modules drawn from modern immunology: somatic hypermutation, idiotypic network regulation, germinal center selection, and microenvironment-driven adaptation.

Performance-critical operations (affinity/distance matrices, clonal selection, network suppression, mutation) are implemented in C++ via [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo), with parallelization support through [BiocParallel](https://bioconductor.org/packages/BiocParallel/).

### Key Features

- **Two tasks** -- clustering and classification on numeric matrices
- **C++ backend** -- BLAS-optimized bulk affinity/distance computation
- **Two APIs** -- functional (`bHIVE()`) for quick use and R6 (`AINet$new()`) for full module composition
- **Multilayer architecture** -- `honeycombHIVE()` for hierarchical prototype refinement across layers
- **Hyperparameter tuning** -- `swarmbHIVE()` with grid search and BiocParallel
- **Gradient refinement** -- `refineB()` post-processing with 5 optimizers and several classification-aware loss functions
- **Composable immune modules** -- mix and match biological mechanisms via dependency injection
- **caret compatible** -- `bHIVEmodel` and `honeycombHIVEmodel` for cross-validation workflows

## Installation

```r
devtools::install_github("BorchLab/bHIVE")
```

## Quick Start

### Functional API

The simplest way to use bHIVE. Works like any R modeling function:

```r
library(bHIVE)
data(iris)
X <- as.matrix(iris[, 1:4])

# Clustering
res <- bHIVE(X, task = "clustering", nAntibodies = 30, maxIter = 20)
table(res$assignments)

# Classification
res <- bHIVE(X, y = iris$Species, task = "classification",
             nAntibodies = 30, maxIter = 20)
table(Predicted = res$assignments, Actual = iris$Species)
```

### R6 API with Modules

For full control, compose an `AINet` with any combination of immune modules:

```r
# Adaptive mutation + idiotypic network regulation
model <- AINet$new(
  nAntibodies = 20,
  maxIter = 30,
  shm = SHMEngine$new(method = "adaptive", base_rate = 0.1),
  idiotypic = IdiotypicNetwork$new(theta_low = 0.01, theta_high = 0.5),
  verbose = FALSE
)
model$fit(X, iris$Species, task = "classification")
table(model$result$assignments)

# Predict on new data
preds <- model$predict(X[1:10, ])
```

## Architecture

### Algorithm

bHIVE evolves a population of antibody vectors to represent structure in data:

1. **Initialization** -- sample from data, random generation, kmeans++, or V(D)J combinatorial assembly
2. **Affinity computation** -- bulk n x m matrix via BLAS (Gaussian, Laplace, polynomial, cosine, Hamming)
3. **Clonal selection + mutation** -- top-k antibodies cloned; mutants generated via configurable SHM strategy
4. **Network regulation** -- suppress redundant antibodies via distance threshold or idiotypic network dynamics
5. **Final assignment** -- data points assigned to nearest antibody by affinity or distance

<img align="center" src="https://github.com/BorchLab/bHIVE/blob/main/www/iterativeGraphic.png">

### Immune Modules

Each module is an R6 class that can be injected into `AINet` via its constructor. All modules are optional -- use only what you need.

| Module | Biological Basis | What It Does |
|:-------|:-----------------|:-------------|
| `SHMEngine` | Somatic hypermutation | 5 mutation strategies: uniform, airs, hotspot, energy, adaptive |
| `IdiotypicNetwork` | Ab-Ab network regulation | Bell-shaped activation dynamics replacing epsilon-threshold suppression |
| `GerminalCenter` | Tfh-B cell interaction | Task-aware quality selection with resource competition |
| `Microenvironment` | Tissue microenvironment cues | Density-dependent zone classification and mutation rate modulation |
| `VDJLibrary` | V(D)J recombination | Combinatorial gene library initialization (PCA, cluster, random partition) |
| `ActivationGate` | Two-signal activation | Costimulatory filtering (density, danger signal, or label entropy) |
| `MemoryPool` | Immunological memory | Archive high-affinity antibodies and recall on distribution shift |
| `ClassSwitcher` | Isotype class switching | IgM (broad) / IgG (specific) / IgA (boundary) kernel width modulation |
| `ConvergentSelector` | Public clonotypes | Cross-repertoire consensus for ensemble methods |

### Multilayer & Tuning

```r
# honeycombHIVE: hierarchical prototype refinement
res <- honeycombHIVE(X, y = iris$Species, task = "classification",
                     layers = 3, nAntibodies = 30,
                     refine = TRUE, refineOptimizer = "adam")

# swarmbHIVE: hyperparameter grid search (parallelizable)
grid <- expand.grid(nAntibodies = c(15, 30), beta = c(3, 5), epsilon = c(0.01, 0.1))
best <- swarmbHIVE(X, y = iris$Species, task = "classification",
                   grid = grid, metric = "accuracy", maxIter = 20)
best$best_params
```

### Gradient Refinement

Fine-tune antibody positions after training with `refineB()`:

```r
res <- bHIVE(X, y = iris$Species, task = "classification",
             nAntibodies = 20, maxIter = 20)

# Adam-based refinement with cross-entropy loss
A_refined <- refineB(res$antibodies, X, y = iris$Species,
                     assignments = res$assignments,
                     task = "classification",
                     loss = "categorical_crossentropy",
                     optimizer = "adam", steps = 10, lr = 0.01)
```

## Affinity & Distance Functions

| Affinity | Formula | Use Case |
|:---------|:--------|:---------|
| `gaussian` | exp(-alpha \|\|x - a\|\|^2) | General purpose (default) |
| `laplace` | exp(-alpha \|\|x - a\|\|) | Heavier tails than Gaussian |
| `polynomial` | (x . a + c)^p | Non-Euclidean similarity |
| `cosine` | (x . a) / (\|\|x\|\| \|\|a\|\|) | Direction-based similarity |
| `hamming` | 1 - (mismatches / d) | Categorical/binary features |

| Distance | Notes |
|:---------|:------|
| `euclidean` | Default; L2 norm |
| `manhattan` | L1 norm |
| `minkowski` | Generalized Lp (parameter p) |
| `cosine` | 1 - cosine similarity |
| `mahalanobis` | Accounts for feature covariance (requires Sigma) |
| `hamming` | Count of differing features |

## Bug Reports / Feature Requests

If you run into any issues or bugs please submit a [GitHub issue](https://github.com/BorchLab/bHIVE/issues) with details of the issue. If possible please include a [reproducible example](https://reprex.tidyverse.org/). Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/BorchLab/bHIVE/issues).

## Contributing

We welcome contributions to the bHIVE project! To contribute:

* Fork the repository.
* Create a feature branch (`git checkout -b feature-branch`).
* Commit your changes (`git commit -m "Add new feature"`).
* Push to the branch (`git push origin feature-branch`).
* Open a pull request.
