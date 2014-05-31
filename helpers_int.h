#ifndef HELPERS_INT_H
#define HELPERS_INT_H
#include <stdbool.h>
#include <stdio.h>

#include "common.h"

#define LIKELY(x) __builtin_expect(!!(x), 1)

#define even(n) ((n) % 2lu == 0)
#define odd(n) ((n) % 2lu == 1)

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
unsigned long int_sqrt(const unsigned long n);

/*
 * Calculate the greatest common divisor of a and b using the Euclidean
 * algorithm.
 */
unsigned long gcd(unsigned long a, unsigned long b);

/*
 * Use the integer square root function to figure out whether a number is a
 * perfect square.
 */
bool is_square(const unsigned long n);

/* Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd. */
void split_int(unsigned long *s, unsigned long *d, const unsigned long n);

/* Compute the jacobi symbol (x/y) of x over y. */
int jacobi(unsigned long x, unsigned long y);

/* This function generates a random integer between low and high. Both low and
 * high are included. */
unsigned long get_random_int(unsigned long low, unsigned long high);

/* Initialises the random number generator using a (static) seed. */
void init_int(void);

#endif
