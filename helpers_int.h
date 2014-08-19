#ifndef HELPERS_INT_H
#define HELPERS_INT_H
#include <stdbool.h>
#include <stdio.h>

#include "common.h"

#define POLY(f) f##_x, f##_1

/* Initialises the random number generator using a (static) seed. */
void init_int(void);

// Replace some GMP functions.
#define even(n) ((n) % 2lu == 0)
#define odd(n) ((n) % 2lu == 1)

/*
 * Calculate the greatest common divisor of a and b using the Euclidean
 * algorithm.
 */
uint64_t gcd(uint64_t a, uint64_t b);

/* Compute the jacobi symbol (x/y) of x over y. */
int jacobi(uint64_t x, uint64_t y);


/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
uint64_t int_sqrt(const uint64_t n);

/*
 * Use the integer square root function to figure out whether a number is a
 * perfect square.
 */
bool is_square(const uint64_t n);


/* This function generates a random integer between low and high. Both low and
 * high are included. */
uint64_t get_random_int(uint64_t low, uint64_t high);


/* Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd. */
void split_int(uint64_t *s, uint64_t *d, const uint64_t n);

#endif
