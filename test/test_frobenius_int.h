/* test/test_frobenius_int.h -- tests for the long long implementation of the QFT
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

#ifndef FROBENIUS_INT_TEST_H
#define FROBENIUS_INT_TEST_H
#include <CUnit/CUnit.h>

void frob_int_int_sqrt(void);
void frob_int_gcd(void);
void frob_int_is_square(void);
void frob_int_jacobi(void);
void frob_int_split(void);
void frob_int_get_random_int(void);
void frob_mult_mod_int(void);
void frob_powm_mod_int(void);
void frob_squares_int(void);
void frob_trial_division_int(void);
void frob_int_rqft_small_primes(void);
void frob_int_rqft_small_composites(void);
void frob_problematic_primes_int(void);
void frob_primelist_int(void);
void frob_larger_primes_int(void);

#endif
