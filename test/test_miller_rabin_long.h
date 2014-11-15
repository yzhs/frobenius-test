/* test/test_miller_rabin_long.h -- tests for the Miller-Rabin test (GMP)
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

#ifndef MILLER_RABIN_LONG_TEST_H
#define MILLER_RABIN_LONG_TEST_H

void mr_some_numbers(void);
void mr_primes(void);
void mr_composites(void);
void mr_composites2(void);
void mr_both(void);

#endif
