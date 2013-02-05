#include "tm.h"

inline double kappa(double *pi, int *population, int pop_size,
                    double mu, R_len_t i, R_len_t j)
{
    return (population[j] - (i == j ? 1 : 0) + mu * pi[j]) /
           (pop_size - 1 + mu);
}

void update_offspring_probs(R_len_t types, int *population, int pop_size,
                            double *offspring_probs)
{
    R_len_t i;

    for (i = 0; i < types; i++) {
        offspring_probs[i] = (double) population[i] / pop_size;
    }
}

void update_ancestor_probs(R_len_t types, double *transitions, double *pi,
                           int *population, int pop_size, double mu,
                           R_len_t offspring_type, double *ancestor_probs)
{
    double p, cumsum;
    R_len_t i;

    /* Compute probability of coalescent event, if allowed */
    ancestor_probs[0] = pop_size > 2 ? population[offspring_type] - 1 : 0;
    cumsum = ancestor_probs[0];

    /* Compute probabilities of mutation events */
    for (i = 0; i < types; i++) {
        p = mu * transitions[offspring_type + types*i] *
            kappa(pi, population, pop_size, mu, offspring_type, i);
        ancestor_probs[i+1] = p;
        cumsum += p;
    }

    /* Normalize vector */
    for (i = 0; i < types + 1; i++) {
        ancestor_probs[i] /= cumsum;
    }
}

void simulate(well1024 *prng, R_len_t types, double *transitions, double *pi,
              int *population, int n, double mu, double *loglik,
              double *correction, int *coalescent_events, int *final_population,
              double *sim_time, double *offspring_probs, double *ancestor_probs)
{
    int pop_size, sen, new_pop_size, subpop_size;
    double start_time, subpop_pi;
    R_len_t offspring_type, ancestor_type, i, coalescent_i;

    pop_size = 0;
    for (i = 0; i < types; i++) {
        final_population[i] = population[i];
        pop_size += population[i];
    }

    start_time = omp_get_wtime();

    *loglik = 0.0;
    sen = 0;
    coalescent_i = 0;
    while (pop_size > n) {
        /* Sample offspring type */
        update_offspring_probs(types, final_population, pop_size,
                               offspring_probs);
        offspring_type = sample(prng, offspring_probs);

        /* Sample ancestor type */
        update_ancestor_probs(types, transitions, pi, final_population,
                              pop_size, mu, offspring_type, ancestor_probs);
        ancestor_type = sample(prng, ancestor_probs) - 1;

        /* Update log-likelihood */
        new_pop_size = pop_size - (ancestor_type < 0 ? 1 : 0);
        *loglik += log(pop_size) + log(pop_size - 1 + mu);
        *loglik -= log(new_pop_size) + log(new_pop_size - 1 + mu);
        if (ancestor_type >= 0) {
            /* Mutation */
            *loglik += log(kappa(pi, final_population, pop_size, mu,
                                 offspring_type, offspring_type));
            *loglik -= log(kappa(pi, final_population, pop_size, mu,
                                 offspring_type, ancestor_type));
            *loglik += log(final_population[ancestor_type] +
                           (ancestor_type != offspring_type ? 1 : 0));
            *loglik -= log(pop_size);
        } else {
            /* Coalescent event */
            *loglik -= log(kappa(pi, final_population, pop_size, mu,
                                 offspring_type, offspring_type));
            *loglik += log(new_pop_size - 1);
            *loglik -= log(final_population[offspring_type]);
        }

        /* Update population size */
        final_population[offspring_type]--;
        if (ancestor_type >= 0) {
            /* Mutation */
            final_population[ancestor_type]++;

            /* Check for the special case of going up to the most recent common
               ancestor, in which case we stop when the population consists of
               two individuals of the same type */
            if (n == 1 && pop_size == 2) {
                break;
            }
        } else {
            /* Coalescent event */
            pop_size--;

            /* Record SEN */
            coalescent_events[coalescent_i++] = sen;
        }

        /* Update SEN */
        sen++;

        R_CheckUserInterrupt();
    }

    /* Compute correction term, if needed */
    *correction = 0.0;
    start_time = omp_get_wtime();
    if (n > 1) {
        *correction += lgammafn(pop_size + 1);
        *correction += lgammafn(mu);
        *correction -= lgammafn(mu + pop_size);
        for (i = 0; i < types; i++) {
            subpop_size = final_population[i];
            subpop_pi = pi[i];
            *correction += lgammafn(subpop_size + mu * subpop_pi);
            *correction -= lgammafn(subpop_size + 1);
            *correction -= lgammafn(mu * subpop_pi);
        }

        *loglik += *correction;
    }

    *sim_time = omp_get_wtime() - start_time;
}

SEXP compute_full_transitions(SEXP unitary_transitions, SEXP loci)
{
    int locus, i_locus, j_locus;
    double first_row_sum, second_row_sum, p;
    R_len_t types, i, j;
    SEXP transitions;

    if (!isMatrix(unitary_transitions)) {
        error("Unitary transitions should be specified as a matrix");
    }
    if (nrows(unitary_transitions) != 2 || ncols(unitary_transitions) != 2) {
        error("The unitary transitions matrix should be 2x2");
    }
    first_row_sum = REAL(unitary_transitions)[0] + REAL(unitary_transitions)[2];
    second_row_sum = REAL(unitary_transitions)[1] + REAL(unitary_transitions)[3];
    if (abs(first_row_sum - 1) > DOUBLE_EPS ||
        abs(second_row_sum - 1) > DOUBLE_EPS) {
        error("Rows of the unitary transitions matrix should sum to 1");
    }
    if (!isInteger(loci) || *INTEGER(loci) <= 0) {
        error("The number of loci should be a positive integer");
    }

    types = 1 << *INTEGER(loci);

    PROTECT(transitions = allocMatrix(REALSXP, types, types));

    for (i = 0; i < types; i++) {
        for (j = 0; j < types; j++) {
            p = 1.0;
            for (locus = 1; locus < types; locus <<= 1) {
                i_locus = (i & locus) == locus;
                j_locus = (j & locus) == locus;
                p *= REAL(unitary_transitions)[i_locus + 2*j_locus];
            }
            REAL(transitions)[j + types*i] = p;
        }
        R_CheckUserInterrupt();
    }

    UNPROTECT(1);

    return transitions;
}

SEXP compute_stationary_distribution(SEXP transitions) {
    double l1norm, l2norm, err;
    R_len_t n, i, j;
    SEXP pi, pi_prev;

    if (!isMatrix(transitions)) {
        error("Transitions should be specified as a matrix");
    }
    if (nrows(transitions) != ncols(transitions)) {
        error("The transition matrix should be square");
    }

    n = nrows(transitions);

    PROTECT(pi = allocVector(REALSXP, n));
    PROTECT(pi_prev = allocVector(REALSXP, n));

    for (i = 0; i < n; i++) {
        REAL(pi)[i] = 1.0 / sqrt(n);
    }

    err = 1.0;
    while (err > DOUBLE_EPS) {
        /* Copy current guess pi into pi_prev */
        for (i = 0; i < n; i++) {
            REAL(pi_prev)[i] = REAL(pi)[i];
        }

        /* Compute new guess pi = transitions pi_prev */
        l1norm = 0.0;
        l2norm = 0.0;
        for (i = 0; i < n; i++) {
            REAL(pi)[i] = 0.0;
            for (j = 0; j < n; j++) {
                REAL(pi)[i] += REAL(transitions)[i + j*n] * REAL(pi_prev)[j];
            }
            l1norm += fabs(REAL(pi)[i]);
            l2norm += R_pow_di(REAL(pi)[i], 2);
        }
        l2norm = sqrt(l2norm);

        /* Normalize pi and update err = (pi_prev-pi)' (pi_prev-pi) */
        err = 0.0;
        for (i = 0; i < n; i++) {
            REAL(pi)[i] /= l2norm;
            err += R_pow_di(REAL(pi_prev)[i] - REAL(pi)[i], 2);
        }

        R_CheckUserInterrupt();
    }

    /* Normalize pi to unit sum */
    for (i = 0; i < n; i++) {
        REAL(pi)[i] /= l1norm;
    }

    UNPROTECT(2);

    return(pi);
}

SEXP estimate_loglik(SEXP transitions_in, SEXP pi_in, SEXP population_in,
                     SEXP n_in, SEXP mu_in, SEXP samples_in, SEXP threads_in)
{
    int *population, n, threads, n_coalescent, *coalescent_events, *final_pops;
    unsigned int *prng_seeds;
    double *transitions, *pi, mu, *offspring_probs, *ancestor_probs, *logliks,
           *corrections, *sim_times, start_time, total_time;
    R_len_t samples, types, i;
    well1024 prng;
    SEXP logliks_vec, corrections_vec, coalescent_events_mat, final_pops_mat,
         sim_times_vec, result, result_names;

    if (!isMatrix(transitions_in) ||
        nrows(transitions_in) != ncols(transitions_in)) {
        error("Transitions should be specified as a square matrix");
    }
    if (!isVector(pi_in)) {
        error("The stationary distribution should be specified as a vector");
    }
    if (length(pi_in) != nrows(transitions_in)) {
        error("The stationary distribution vector should have as many entries as there are types");
    }
    if (!isVector(population_in)) {
        error("The initial population should be specified as a vector");
    }
    if (length(population_in) != nrows(transitions_in)) {
        error("The initial population vector should have as many entries as there are types");
    }

    types = nrows(transitions_in);
    transitions = REAL(transitions_in);
    pi = REAL(pi_in);
    population = INTEGER(population_in);
    n = asInteger(n_in);
    mu = asReal(mu_in);
    samples = asInteger(samples_in);

    if (n <= 0) {
        error("The target population size should be a positive integer");
    }
    if (mu <= 0) {
        error("The mutation rate should be positive");
    }

    /* Set number of threads, if requested */
    if (!isNull(threads_in)) {
        omp_set_num_threads(asInteger(threads_in));
    }

    /* Compute number of coalescent events */
    n_coalescent = -n;
    for (i = 0; i < types; i++) {
        n_coalescent += population[i];
    }
    if (n == 1) {
        n_coalescent--;
    }

    PROTECT(logliks_vec = allocVector(REALSXP, samples));
    PROTECT(corrections_vec = allocVector(REALSXP, samples));
    PROTECT(coalescent_events_mat = allocMatrix(INTSXP, n_coalescent, samples));
    PROTECT(final_pops_mat = allocMatrix(INTSXP, types, samples));
    PROTECT(sim_times_vec = allocVector(REALSXP, samples));

    logliks = REAL(logliks_vec);
    corrections = REAL(corrections_vec);
    coalescent_events = INTEGER(coalescent_events_mat);
    final_pops = INTEGER(final_pops_mat);
    sim_times = REAL(sim_times_vec);

    #pragma omp parallel \
     shared (transitions, pi, population, n, threads, mu, samples, logliks, \
             corrections, coalescent_events, final_pops, sim_times, \
             start_time, total_time, prng_seeds) \
     private (i, prng, offspring_probs, ancestor_probs)
    {
        #pragma omp single
        {
            well1024_init(&prng, time_seed());
            threads = omp_get_num_threads();
            prng_seeds = (unsigned int *) malloc(sizeof(unsigned int) * threads);
            for (i = 0; i < threads; i++) {
                prng_seeds[i] = well1024_uint_rand(&prng);
            }
        }

        well1024_init(&prng, prng_seeds[omp_get_thread_num()]);
        offspring_probs = (double *) malloc(sizeof(double) * types);
        ancestor_probs = (double *) malloc(sizeof(double) * (types + 1));

        #pragma omp single
        {
            start_time = omp_get_wtime();
        }

        #pragma omp for
        for (i = 0; i < samples; i++) {
            simulate(&prng, types, transitions, pi, population, n, mu,
                     &logliks[i], &corrections[i],
                     &coalescent_events[i*n_coalescent],
                     &final_pops[i*types], &sim_times[i], offspring_probs,
                     ancestor_probs);
        }

        #pragma omp single
        {
            total_time = omp_get_wtime() - start_time;
        }

        free(offspring_probs);
        free(ancestor_probs);
    }

    free(prng_seeds);

    PROTECT(result = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(result, 0, population_in);
    SET_VECTOR_ELT(result, 1, n_in);
    SET_VECTOR_ELT(result, 2, mu_in);
    SET_VECTOR_ELT(result, 3, logliks_vec);
    SET_VECTOR_ELT(result, 4, corrections_vec);
    SET_VECTOR_ELT(result, 5, coalescent_events_mat);
    SET_VECTOR_ELT(result, 6, final_pops_mat);
    SET_VECTOR_ELT(result, 7, sim_times_vec);
    SET_VECTOR_ELT(result, 8, ScalarReal(total_time));

    PROTECT(result_names = allocVector(VECSXP, 9));
    SET_VECTOR_ELT(result_names, 0, mkChar("population"));
    SET_VECTOR_ELT(result_names, 1, mkChar("n"));
    SET_VECTOR_ELT(result_names, 2, mkChar("mu"));
    SET_VECTOR_ELT(result_names, 3, mkChar("logliks"));
    SET_VECTOR_ELT(result_names, 4, mkChar("corrections"));
    SET_VECTOR_ELT(result_names, 5, mkChar("coalescent.events"));
    SET_VECTOR_ELT(result_names, 6, mkChar("final.populations"));
    SET_VECTOR_ELT(result_names, 7, mkChar("simulation.times"));
    SET_VECTOR_ELT(result_names, 8, mkChar("total.time"));
    setAttrib(result, R_NamesSymbol, result_names);

    UNPROTECT(7);

    return result;
}

