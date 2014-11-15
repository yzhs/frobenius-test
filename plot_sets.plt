# plot_mults.plt -- plot runtime per input set
#
# Copyright 2014 by Colin Benner <colin-software@yzhs.de>
#
# This file is part of frobenius-test.
#
# frobenius-test is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# frobenius-test is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.
set datafile separator ','

set terminal epslatex color size 17cm,10cm

set style line 1 lt rgb "red" lw 2
set style line 2 lt rgb "orange" lw 2
set style line 3 lt rgb "#f0f000" lw 2
set style line 4 lt rgb "green" lw 2
set style line 5 lt rgb "cyan" lw 2
set style line 6 lt rgb "blue" lw 2
set style line 7 lt rgb "violet" lw 2

set xlabel '$\lfloor\log_2 n\rfloor$'
set ylabel 'Laufzeit in Sekunden'

set format x '$2^{%L}$'

set xrange [32:131072]
set logscale x 2
set logscale y

set key top left
set ytics 1e-6,1e2,1e8

# Full test
set output 'pic/primes.tex'
plot 'data/primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/primes_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/primes_frob.csv' using 1:2 title 'Frobenius'

set output 'pic/all_primes.tex'
plot 'data/primes_gmp.csv'           using 1:2 ls 1 title 'GMP ($2^k+\varepsilon$)', \
     'data/primes_mr.csv'            using 1:2 ls 2 title 'Miller-Rabin ($2^k+\varepsilon$)', \
     'data/primes_frob.csv'          using 1:2 ls 3 title 'Frobenius ($2^k+\varepsilon$)', \
     'data/mersenne_primes_gmp.csv'  using 1:2 ls 4 title 'GMP ($2^p-1$)', \
     'data/mersenne_primes_mr.csv'   using 1:2 ls 5 title 'Miller-Rabin ($2^p-1$)', \
     'data/mersenne_primes_frob.csv' using 1:2 ls 6 title 'Frobenius ($2^p-1$)'

set output 'pic/all_primes_errorbars.tex'
plot 'data/primes_gmp.csv'           using 1:2:3:4 with errorbars ls 1 title 'GMP ($2^k+\varepsilon$)', \
     'data/primes_mr.csv'            using 1:2:3:4 with errorbars ls 2 title 'Miller-Rabin ($2^k+\varepsilon$)', \
     'data/primes_frob.csv'          using 1:2:3:4 with errorbars ls 3 title 'Frobenius ($2^k+\varepsilon$)', \
     'data/mersenne_primes_gmp.csv'  using 1:2:3:4 with errorbars ls 4 title 'GMP ($2^p-1$)', \
     'data/mersenne_primes_mr.csv'   using 1:2:3:4 with errorbars ls 5 title 'Miller-Rabin ($2^p-1$)', \
     'data/mersenne_primes_frob.csv' using 1:2:3:4 with errorbars ls 6 title 'Frobenius ($2^p-1$)'


set output 'pic/composites.tex'
plot 'data/composites_gmp.csv'  using 1:2 title 'GMP', \
     'data/composites_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/composites_frob.csv' using 1:2 title 'Frobenius'


set output 'pic/mersenne_primes.tex'
plot 'data/mersenne_primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/mersenne_primes_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/mersenne_primes_frob.csv' using 1:2 title 'Frobenius'

set xrange [32:2048]
set output 'pic/mersenne_numbers.tex'
plot 'data/mersenne_numbers_gmp.csv'  using 1:2 title 'GMP', \
     'data/mersenne_numbers_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/mersenne_numbers_frob.csv' using 1:2 title 'Frobenius'
set xrange [32:131072]


# Precomputation
set output 'pic/prep_primes.tex'
plot 'data/prep_primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/prep_primes_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/prep_primes_frob.csv' using 1:2 title 'Frobenius'


set output 'pic/prep_composites.tex'
plot 'data/prep_composites_gmp.csv'  using 1:2 title 'GMP', \
     'data/prep_composites_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/prep_composites_frob.csv' using 1:2 title 'Frobenius'


set output 'pic/prep_mersenne_primes.tex'
plot 'data/prep_mersenne_primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/prep_mersenne_primes_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/prep_mersenne_primes_frob.csv' using 1:2 title 'Frobenius'

set xrange [32:2048]
set output 'pic/prep_mersenne_numbers.tex'
plot 'data/prep_mersenne_numbers_gmp.csv'  using 1:2 title 'GMP', \
     'data/prep_mersenne_numbers_mr.csv'   using 1:2 title 'Miller-Rabin', \
     'data/prep_mersenne_numbers_frob.csv' using 1:2 title 'Frobenius'
set xrange [32:131072]


set ylabel "Laufzeit relativ zu Miller-Rabin"
set grid y

set key top right
set yrange [1:80]
set ytics 1,2,80

# Normalized runtime
set output 'pic/all_primes_normalized.tex'
plot 'data/normalized_primes_gmp.csv'           using 1:2 title 'GMP ($2^k+\varepsilon$)', \
     'data/normalized_primes_frob.csv'          using 1:2 title 'Frobenius ($2^k+\varepsilon$)', \
     'data/normalized_mersenne_primes_gmp.csv'  using 1:2 title 'GMP ($2^p-1$)', \
     'data/normalized_mersenne_primes_frob.csv' using 1:2 title 'Frobenius ($2^p-1$)'

set output 'pic/normalized_primes.tex'
plot 'data/normalized_primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/normalized_primes_frob.csv' using 1:2 title 'Frobenius'


set output 'pic/normalized_composites.tex'
plot 'data/normalized_composites_gmp.csv'  using 1:2 title 'GMP', \
     'data/normalized_composites_frob.csv' using 1:2 title 'Frobenius'


set output 'pic/normalized_mersenne_primes.tex'
plot 'data/normalized_mersenne_primes_gmp.csv'  using 1:2 title 'GMP', \
     'data/normalized_mersenne_primes_frob.csv' using 1:2 title 'Frobenius'

set xrange [32:2048]
set output 'pic/normalized_mersenne_numbers.tex'
plot 'data/normalized_mersenne_numbers_gmp.csv'  using 1:2 title 'GMP', \
     'data/normalized_mersenne_numbers_frob.csv' using 1:2 title 'Frobenius'
set xrange [32:131072]
