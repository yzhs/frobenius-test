#ifndef HELPERS_INT_H
#define HELPERS_INT_H
#include <stdbool.h>

#define LIKELY(x) __builtin_expect(!!(x), 1)

#define die(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

#ifdef DEBUG
#define debug(...) printf(__VA_ARGS__)
#else
#define debug(...)
#endif

#define len(a) (sizeof(a)/sizeof(a[0]))

#define even(n) ((n) % 2 == 0)
#define odd(n) ((n) % 2 == 1)

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
unsigned long int_sqrt(unsigned long n);

/*
 * Calculate the greatest common divisor of a and b using the Euclidean
 * algorithm.
 */
unsigned long gcd(unsigned long a, unsigned long b);

/*
 * Use the integer square root function to figure out whether a number is a
 * perfect square.
 */
bool is_square(unsigned long n);

/* Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd. */
unsigned long split(unsigned long *d, unsigned long n);

/*
 * This function generates a random integer between 2 and n-2.  This function
 * will fail if n is equal to 3 and produce weird results for smaller ns.
 */
unsigned long get_random(unsigned long n);

/* Compute the jacobi symbol (x/y) of x over y. */
int jacobi(unsigned long x, unsigned long y);

/* This function generates a random integer between 2 and n-2. */
unsigned long get_random(unsigned long n);

/* Initialises the random number generator using a (static) seed. */
void init();

#endif
