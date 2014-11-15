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
set output 'plots/primes.tex'
plot 'processed/primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/primes_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/primes_frob.csv' using 1:2 ls 3 title 'Frobenius'

set output 'plots/all_primes.tex'
plot 'processed/primes_gmp.csv'           using 1:2 ls 1 title 'GMP ($2^k+\varepsilon$)', \
     'processed/primes_mr.csv'            using 1:2 ls 2 title 'Miller-Rabin ($2^k+\varepsilon$)', \
     'processed/primes_frob.csv'          using 1:2 ls 3 title 'Frobenius ($2^k+\varepsilon$)', \
     'processed/mersenne_primes_gmp.csv'  using 1:2 ls 4 title 'GMP ($2^p-1$)', \
     'processed/mersenne_primes_mr.csv'   using 1:2 ls 5 title 'Miller-Rabin ($2^p-1$)', \
     'processed/mersenne_primes_frob.csv' using 1:2 ls 6 title 'Frobenius ($2^p-1$)'

set output 'plots/all_primes_errorbars.tex'
plot 'processed/primes_gmp.csv'           using 1:2:3:4 with errorbars ls 1 title 'GMP ($2^k+\varepsilon$)', \
     'processed/primes_mr.csv'            using 1:2:3:4 with errorbars ls 2 title 'Miller-Rabin ($2^k+\varepsilon$)', \
     'processed/primes_frob.csv'          using 1:2:3:4 with errorbars ls 3 title 'Frobenius ($2^k+\varepsilon$)', \
     'processed/mersenne_primes_gmp.csv'  using 1:2:3:4 with errorbars ls 4 title 'GMP ($2^p-1$)', \
     'processed/mersenne_primes_mr.csv'   using 1:2:3:4 with errorbars ls 5 title 'Miller-Rabin ($2^p-1$)', \
     'processed/mersenne_primes_frob.csv' using 1:2:3:4 with errorbars ls 6 title 'Frobenius ($2^p-1$)'


set output 'plots/composites.tex'
plot 'processed/composites_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/composites_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/composites_frob.csv' using 1:2 ls 3 title 'Frobenius'


set output 'plots/mersenne_primes.tex'
plot 'processed/mersenne_primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/mersenne_primes_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/mersenne_primes_frob.csv' using 1:2 ls 3 title 'Frobenius'

set xrange [32:2048]
set output 'plots/mersenne_numbers.tex'
plot 'processed/mersenne_numbers_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/mersenne_numbers_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/mersenne_numbers_frob.csv' using 1:2 ls 3 title 'Frobenius'
set xrange [32:131072]


# Precomputation
set output 'plots/prep_primes.tex'
plot 'processed/prep_primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/prep_primes_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/prep_primes_frob.csv' using 1:2 ls 3 title 'Frobenius'


set output 'plots/prep_composites.tex'
plot 'processed/prep_composites_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/prep_composites_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/prep_composites_frob.csv' using 1:2 ls 3 title 'Frobenius'


set output 'plots/prep_mersenne_primes.tex'
plot 'processed/prep_mersenne_primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/prep_mersenne_primes_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/prep_mersenne_primes_frob.csv' using 1:2 ls 3 title 'Frobenius'

set xrange [32:2048]
set output 'plots/prep_mersenne_numbers.tex'
plot 'processed/prep_mersenne_numbers_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/prep_mersenne_numbers_mr.csv'   using 1:2 ls 2 title 'Miller-Rabin', \
     'processed/prep_mersenne_numbers_frob.csv' using 1:2 ls 3 title 'Frobenius'
set xrange [32:131072]


set ylabel "Laufzeit relativ zu Miller-Rabin"
set grid y

set key top right
set yrange [1:80]
set ytics 1,2,80

# Normalized runtime
set output 'plots/all_primes_normalized.tex'
plot 'processed/normalized_primes_gmp.csv'           using 1:2 ls 1 title 'GMP ($2^k+\varepsilon$)', \
     'processed/normalized_primes_frob.csv'          using 1:2 ls 2 title 'Frobenius ($2^k+\varepsilon$)', \
     'processed/normalized_mersenne_primes_gmp.csv'  using 1:2 ls 3 title 'GMP ($2^p-1$)', \
     'processed/normalized_mersenne_primes_frob.csv' using 1:2 ls 4 title 'Frobenius ($2^p-1$)'

set output 'plots/normalized_primes.tex'
plot 'processed/normalized_primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/normalized_primes_frob.csv' using 1:2 ls 2 title 'Frobenius'


set output 'plots/normalized_composites.tex'
plot 'processed/normalized_composites_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/normalized_composites_frob.csv' using 1:2 ls 2 title 'Frobenius'


set output 'plots/normalized_mersenne_primes.tex'
plot 'processed/normalized_mersenne_primes_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/normalized_mersenne_primes_frob.csv' using 1:2 ls 2 title 'Frobenius'

set xrange [32:2048]
set output 'plots/normalized_mersenne_numbers.tex'
plot 'processed/normalized_mersenne_numbers_gmp.csv'  using 1:2 ls 1 title 'GMP', \
     'processed/normalized_mersenne_numbers_frob.csv' using 1:2 ls 2 title 'Frobenius'
set xrange [32:131072]

# vim: set ft=gnuplot :
