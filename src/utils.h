#ifndef UTILS_H
#define UTILS_H

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <Rinternals.h>

#define WELL_MAT3POS(t, v) (v^(v>>t))
#define WELL_MAT3NEG(t, v) (v^(v<<(-t)))

typedef struct {
    unsigned int seed;
    unsigned state_n;
    unsigned int state[32];
} well1024;

unsigned int time_seed();
void well1024_init(well1024 *, unsigned int);
double well1024_unif_rand(well1024 *);
unsigned int well1024_uint_rand(well1024 *);
R_len_t sample(well1024 *, double *);

#endif

