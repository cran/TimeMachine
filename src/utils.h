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

