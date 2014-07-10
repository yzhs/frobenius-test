#ifndef SAMLL_PRIMES_H
#define SAMLL_PRIMES_H

// Statically determined list of all primes less than 50000.
#if B > 44958
extern const unsigned short prime_list[5133];
#else
extern const unsigned short prime_list[4670];
#endif

#endif
