#ifndef BHIVE_SHM_H
#define BHIVE_SHM_H

// Somatic Hypermutation strategy helpers shared between shm.cpp (the
// standalone shm_mutate_cpp entry point) and clonal_selection.cpp (which
// dispatches per-clone mutation through `mutate_by_method`).

#include <RcppArmadillo.h>
#include "affinity_distance.h"

arma::rowvec mutate_uniform(const arma::rowvec& antibody, double affinity,
                            int iter, double decay, double mutationMin);

arma::rowvec mutate_airs(const arma::rowvec& antibody, double affinity,
                        double c_rate, double temperature);

arma::rowvec mutate_hotspot(const arma::rowvec& antibody,
                            const arma::rowvec& data_point,
                            double affinity, double base_rate);

arma::rowvec mutate_energy(const arma::rowvec& antibody, double affinity,
                          double E_0);

// Adaptive (Adam-style) mutation. Updates m1_state.row(j) and
// m2_state.row(j) in place. Caller is responsible for ensuring
// m1_state/m2_state have shape (nAntibodies, nFeatures).
arma::rowvec mutate_adaptive(const arma::rowvec& antibody,
                            const arma::rowvec& data_point,
                            double affinity, int iter, double mutationMin,
                            double base_rate, double beta1, double beta2,
                            double adam_epsilon,
                            arma::mat& m1_state, arma::mat& m2_state,
                            arma::uword j);

// Dispatch a single mutation by method name. For methods that do not need
// data_point (uniform, airs, energy), the value is ignored. For methods
// that do not need m1/m2 state, the matrices are ignored.
arma::rowvec mutate_by_method(const std::string& method,
                              const arma::rowvec& antibody,
                              const arma::rowvec& data_point,
                              double affinity, int iter,
                              double decay, double mutationMin,
                              double c_rate, double temperature,
                              double E_0, double base_rate,
                              double beta1, double beta2,
                              double adam_epsilon,
                              arma::mat& m1_state, arma::mat& m2_state,
                              arma::uword j);

#endif
