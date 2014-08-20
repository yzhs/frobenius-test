#ifndef COMMON_H
#define COMMON_H

#include <stdarg.h>
#include <stdint.h>

#define LIKELY(x) __builtin_expect(!!(x), 1)

extern int enable_logging;

//#define die(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

// Print an error message and terminate the programm.
int die(const char *format, ...);

// Print a debugging message if debugging output is enabled.
#define debug(...) do { if (enable_logging) fprintf(stderr, __VA_ARGS__); } while (0)


// Upper bound for the prime numbers to be use for trial division.
#define B 44958lu

// Possible result of primality tests that prove compositeness.
typedef enum { composite = 0, probably_prime, prime } Primality;


// Get the number of elements a (statically sized) array has.
#define len(a) (sizeof(a) / sizeof(a[0]))

// Return a certain value x after deallocating the local big integer variables.
#define ret(x) do { result = (x); goto exit; } while (0)

#define TEST_DATA_PATH "test/data/"

#endif
