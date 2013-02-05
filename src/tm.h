#ifndef TM_H
#define TM_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h>

#include "utils.h"

inline double kappa(double *, int *, int, double, R_len_t, R_len_t);
void update_offspring_probs(R_len_t, int *, int, double *);
void update_ancestor_probs(R_len_t, double *, double *, int *, int, double,
                           R_len_t, double *);
void simulate(well1024 *, R_len_t, double *, double *, int *, int, double,
              double *, double *, int *, int *, double *, double *, double *);

SEXP compute_full_transitions(SEXP, SEXP);
SEXP compute_stationary_distribution(SEXP);
SEXP estimate_loglik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif

