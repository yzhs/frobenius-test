/* test/test_frobenius_long.h -- tests for the GMP implementation of the QFT
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
#ifndef TEST_FROBENIUS_LONG_H
#define TEST_FROBENIUS_LONG_H
#include <CUnit/CUnit.h>

void frob_mult_x(void);
void frob_sigma_basics(void);
void frob_sigma_short_integer(void);
void frob_sigma_power(void);
void frob_power_basics(void);
void frob_power_x_lucas(void);
void frob_inverse(void);
void frob_fast_algorithm1(void);
void frob_fast_algorithm2(void);
void frob_fast_algorithm3(void);
void frob_split(void);
void frob_get_random(void);
void frob_mult_mod(void);
void frob_powm_mod(void);
void frob_squares(void);
void frob_trial_division(void);
void frob_rqft_small_primes(void);
void frob_rqft_small_composites(void);
void frob_problematic_primes(void);
void frob_primelist(void);
void frob_larger_primes(void);
void frob_composites(void);


#endif
