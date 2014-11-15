/* frobenius_int.h -- the quadratic frobenius test (long long version)
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

#ifndef FROBENIUS_INT_H
#define FROBENIUS_INT_H

#include "common.h"

Primality steps_3_4_5_int(const uint64_t n, const uint64_t b, const uint64_t c);
Primality QFT_int(const unsigned n, const unsigned b, const unsigned c);
Primality RQFT_int(const unsigned n, const unsigned k);

#endif
