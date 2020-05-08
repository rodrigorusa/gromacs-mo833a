/*
 * This file is just instrumentation to capture paramount iterations time
 *
 */
#ifndef GMX_UTILITY_TIME_H
#define GMX_UTILITY_TIME_H

#include <sys/time.h>

extern double T_START_MAIN;

// Function to get time in seconds
inline static double mysecond() {
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

#endif
