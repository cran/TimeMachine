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

#include "utils.h"

unsigned int time_seed()
{
    unsigned char *p;
    unsigned int seed;
    size_t i;
    time_t now;

    now = time(NULL);
    p = (unsigned char *) &now;
    seed = 0;
    for (i = 0; i < sizeof(now); i++) {
        seed = seed * (UCHAR_MAX + 2U) + p[i];
    }

    return seed;
}

void well1024_init(well1024 *prng, unsigned int seed)
{
    int i;

    prng->seed = seed;
    prng->state_n = 0;
    prng->state[0] = seed & 0xFFFFFFFFU;
    for (i = 1; i < 32; i++) {
        prng->state[i] = (69069 * prng->state[i-1]) & 0xFFFFFFFFU;
    }
}

double well1024_unif_rand(well1024 *prng)
{
    unsigned int *state, state_n, z0, z1, z2;

    state = prng->state;
    state_n = prng->state_n;
    z0 = state[(state_n + 31) & 0x0000001FUL];
    z1 = state[state_n] ^ WELL_MAT3POS(8, state[(state_n + 3) & 0x0000001FUL]);
    z2 = WELL_MAT3NEG(-19, state[(state_n + 24) & 0x0000001FUL]) ^
         WELL_MAT3NEG(-14, state[(state_n + 10) & 0x0000001FUL]);
    state[state_n] = z1 ^ z2;
    state[(state_n + 31) & 0x0000001FUL] = WELL_MAT3NEG(-11, z0) ^
                                           WELL_MAT3NEG( -7, z1) ^
                                           WELL_MAT3NEG(-13, z2);
    prng->state_n = (state_n + 31) & 0x0000001FUL;

    return ((double) state[prng->state_n] * 2.32830643653869628906E-10);
}

unsigned int well1024_uint_rand(well1024 *prng)
{
    return (unsigned int) (well1024_unif_rand(prng) * UINT_MAX);
}

R_len_t sample(well1024 *prng, double *probabilities)
{
    double u, cumsum;
    R_len_t i;

    u = well1024_unif_rand(prng);

    cumsum = probabilities[0];
    i = 0;
    while (u > cumsum) {
        cumsum += probabilities[++i];
    }
    return i;
}

