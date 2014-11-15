/* frobenius.h -- the quadratic frobenius test (GMP version)
 *
 * Copyright 2014 by Colin Benner <colin-software@yzhs.de>
 *
 * This file is part of frobenius-test.
 *
 * frobenius-test is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * frobenius-test is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gmp.h>

#include "common.h"

Primality steps_1_2(const mpz_t n);
Primality steps_3_4_5(const mpz_t n, const mpz_t b, const mpz_t c);

Primality QFT(const mpz_t n, const mpz_t b, const mpz_t c);
Primality RQFT(const mpz_t n, const unsigned k);

#endif

