/*
 * This file is part of the TimeMachine.
 * Copyright (C) 2013 Gianluca Campanella <gianluca@campanella.org>
 *
 * The TimeMachine is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * The TimeMachine is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * theTimeMachine. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TM_H
#define TM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <R.h>
#include <Rmath.h>

#include "utils.h"

static inline double kappa(double *, int *, int, double, R_len_t, R_len_t);
void update_offspring_probs(R_len_t, int *, int, double *);
void update_ancestor_probs(R_len_t, double *, double *, int *, int, double,
                           R_len_t, double *);
void simulate(well1024 *, R_len_t, double *, double *, int *, int, double,
              double *, double *, int *, int *, double *, double *, double *);

SEXP compute_full_transitions(SEXP, SEXP);
SEXP compute_stationary_distribution(SEXP);
SEXP estimate_loglik(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

#endif

