# bHIVE 0.99.2

## Breaking Changes
* Removed regression task support across the package. The previous regression
  implementation was unreliable and is being redesigned. `bHIVE()`,
  `honeycombHIVE()`, `swarmbHIVE()`, `refineB()`, `AINet`, `GerminalCenter`,
  `bHIVEmodel`, `honeycombHIVEmodel`, and `visualizeHIVE()` no longer accept
  `task = "regression"` or numeric `y`.
* Removed regression-only loss functions (`mse`, `huber`, `poisson`) and
  regression metrics (`rmse`, `mae`, `r2`) from `swarmbHIVE()`.
* Removed `refineHuberDelta` parameter from `honeycombHIVE()` and
  `honeycombHIVEmodel`.
* `final_assignment_cpp()` no longer accepts `antibody_values` or
  `overall_mean` arguments.

# bHIVE 0.99.1

## C++ Backend
* Added RcppArmadillo backend for BLAS-optimized bulk affinity and distance
  matrix computation, replacing per-element R loops
* Scalar affinity function for clone/mutate hot path avoids 1x1 matrix
  allocation overhead
* C++ implementations for clonal selection iteration, network suppression,
  kmeans++ initialization, final assignment, somatic hypermutation (5 methods),
  and idiotypic network dynamics

## R6 Class Architecture
* New `ImmuneRepertoire` class for antibody collections with metadata tracking
  (isotype, state, age, lineage)
* New `ImmuneAlgorithm` abstract base class with fit/predict/summary interface
* New `AINet` class wrapping core algorithm with composable module injection
* Both R6 composition and functional API equally supported

## New Modules
* `SHMEngine` — Five somatic hypermutation strategies: uniform (original
  behavior), airs (affinity-proportional), hotspot (feature-gradient-weighted),
  energy (budget-constrained), and adaptive (per-feature Adam-like moment
  tracking)
* `IdiotypicNetwork` — Antibody-antibody network dynamics with bell-shaped
  activation function replacing epsilon-threshold suppression
* `GerminalCenter` — T-follicular helper mediated selection with task-aware
  quality scoring and resource competition
* `Microenvironment` — Density-dependent zone classification
  (stable/explore/boundary) with chemokine-like gradient computation
* `VDJLibrary` — Combinatorial V(D)J gene library initialization via PCA,
  k-means clustering, or random partition of feature space
* `ActivationGate` — Two-signal activation gate requiring both antigen
  recognition and costimulatory context (density, danger, or entropy)
* `MemoryPool` — Archive high-affinity antibodies as long-lived memory cells
  with threshold-based recall
* `ClassSwitcher` — Isotype class switching (IgM broad, IgG specific, IgA
  boundary) modulating effective kernel width
* `ConvergentSelector` — Cross-repertoire consensus identification of public
  antibodies for ensemble methods

## Documentation
* Complete README rewrite covering functional API, R6 API, module reference
  table, and architecture overview
* New pkgdown website with Bootstrap 5 Flatly theme, organized reference groups,
  and tutorial navigation
* New article: "Composing Immune Modules" — R6 composition patterns for all 9
  modules with worked examples
* New article: "Advanced Tuning & Workflows" — swarmbHIVE grid search,
  honeycombHIVE multilayer refinement, refineB optimizer comparison, caret
  integration, and visualizeHIVE plot types
* New article: "Algorithm & Biological Foundations" — comprehensive mathematical
  reference covering all affinity kernels, distance functions, SHM strategies,
  idiotypic ODE system, germinal center selection, and parameter guidance
* Added roxygen @examples to refineB, bHIVEmodel, and ImmuneAlgorithm (now 90%
  example coverage)
* GitHub Actions workflow for automated pkgdown deployment to borch.dev

## Package Infrastructure
* Created `R/bHIVE-package.R` with roxygen-managed `@useDynLib` and
 `@importFrom Rcpp sourceCpp` directives
* Added `%||%` operator `@name null-coalesce` to avoid illegal characters in Rd
  `\name` field
* Moved tutorial vignettes to `vignettes/articles/` (pkgdown-only, not installed
  with package) to reduce installed size
* Added pkgdown configuration (`_pkgdown.yml`) with 7 reference groups and
  structured article hierarchy

## Testing
* Added comprehensive unit tests for all 12 R6 module classes (ImmuneRepertoire,
  ImmuneAlgorithm, AINet, SHMEngine, IdiotypicNetwork, GerminalCenter,
  Microenvironment, VDJLibrary, ActivationGate, MemoryPool, ClassSwitcher,
  ConvergentSelector) and C++ backend functions
* Test coverage increased from ~26% to ~85% (681 tests, 0 failures)

## BiocCheck Compliance
* Replaced all `sapply()` calls with `vapply()` in bHiVE.R and visualizeHIVE.R
* Replaced all `1:n` patterns with `seq_len()` / `seq_along()`
* Removed `install.packages()` calls from vignettes
* Removed `LazyData: true` from DESCRIPTION
* Updated R dependency to >= 4.5.0
* Updated biocViews to `Software, Clustering, Classification, Regression,
  Network`
* Added class-level `@param` documentation for all R6 initialize() arguments
  across 11 module classes
* Added `@param ... Not used.` to all R6 `print()` methods

## Bug Fixes
* Fixed kmeans++ sampling: corrected cumulative sum fallback index and
  runif-to-integer cast
* Fixed cosine similarity epsilon (1e-12) for numerical stability
* Fixed division-by-zero guard in idiotypic dynamics when theta_low equals
  theta_high
* Fixed VDJLibrary NaN in kmeans for single-allele edge cases
* Fixed VDJLibrary NA subscript when feature dimensions not divisible by 3
* Removed duplicate `Classification` from biocViews, added `Clustering`
* Fixed vignette YAML parsing error from bare `---` horizontal rule

# bHIVE 0.99.0

* Initial submission version with core AIS functionality
* Clonal selection, network suppression, and mutation for clustering,
  classification, and regression
* honeycombHIVE multilayer architecture
* swarmbHIVE hyperparameter tuning via BiocParallel
* caret model integration (bHIVEmodel, honeycombHIVEmodel)
* Visualization utilities via ggplot2
