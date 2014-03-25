#ifndef COMMON_H
#define COMMON_H

#define die(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

#define len(a) (sizeof(a)/sizeof(a[0]))

typedef enum { composite = 0, probably_prime, prime } Primality;

#endif