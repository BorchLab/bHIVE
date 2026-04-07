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

## Bug Fixes
* Fixed kmeans++ sampling: corrected cumulative sum fallback index and
  runif-to-integer cast
* Fixed cosine similarity epsilon (1e-12) for numerical stability
* Fixed division-by-zero guard in idiotypic dynamics when theta_low equals
  theta_high
* Fixed VDJLibrary NaN in kmeans for single-allele edge cases
* Fixed VDJLibrary NA subscript when feature dimensions not divisible by 3
* Removed duplicate `Classification` from biocViews, added `Clustering`

# bHIVE 0.99.0

* Initial submission version with core AIS functionality
* Clonal selection, network suppression, and mutation for clustering,
  classification, and regression
* honeycombHIVE multilayer architecture
* swarmbHIVE hyperparameter tuning via BiocParallel
* caret model integration (bHIVEmodel, honeycombHIVEmodel)
* Visualization utilities via ggplot2
