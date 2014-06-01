#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <gmp.h>

#include "common.h"

Primality miller_rabin(const mpz_t n, const unsigned k);

#endif


