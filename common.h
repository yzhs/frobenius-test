/* common.h -- functionality for all algorithms and versions (long long and GMP)
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

#ifndef COMMON_H
#define COMMON_H

#include <stdarg.h>
#include <stdint.h>

#define LIKELY(x) __builtin_expect(!!(x), 1)

extern int enable_logging;

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
