# plot_algorithms.plt -- plot runtime per algorithm
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
set style line 3 lt rgb "#00B000" lw 2
set style line 4 lt rgb "black" lw 2
set style line 5 lt rgb "cyan" lw 2
set style line 6 lt rgb "blue" lw 2
set style line 7 lt rgb "violet" lw 2

set xlabel '$\log_2 n$'
set ylabel 'Laufzeit in Sekunden'

set format x '$2^{%L}$'

set xrange [32:131072]
set logscale x 2
set logscale y

set key top left
set ytics 1e-6,1e2,1e8

# Full test
set output 'plots/gmp.tex'
plot 'processed/primes_gmp.csv'            using 1:2 ls 1 title 'Primzahlen', \
     'processed/composites_gmp.csv'        using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/mersenne_numbers_gmp.csv'  using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/mersenne_primes_gmp.csv'   using 1:2 ls 4 title 'Mersenne-Primzahlen'

set output 'plots/mr.tex'
plot 'processed/primes_mr.csv'             using 1:2 ls 1 title 'Primzahlen', \
     'processed/composites_mr.csv'         using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/mersenne_numbers_mr.csv'   using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/mersenne_primes_mr.csv'    using 1:2 ls 4 title 'Mersenne-Primzahlen'

set output 'plots/frob.tex'
plot 'processed/primes_frob.csv'           using 1:2 ls 1 title 'Primzahlen', \
     'processed/composites_frob.csv'       using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/mersenne_numbers_frob.csv' using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/mersenne_primes_frob.csv'  using 1:2 ls 4 title 'Mersenne-Primzahlen'


# Just the precomputation
set output 'plots/prep_gmp.tex'
plot 'processed/prep_primes_gmp.csv'            using 1:2 ls 1 title 'Primzahlen', \
     'processed/prep_composites_gmp.csv'        using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/prep_mersenne_numbers_gmp.csv'  using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/prep_mersenne_primes_gmp.csv'   using 1:2 ls 4 title 'Mersenne-Primzahlen'

set output 'plots/prep_mr.tex'
plot 'processed/prep_primes_mr.csv'             using 1:2 ls 1 title 'Primzahlen', \
     'processed/prep_composites_mr.csv'         using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/prep_mersenne_numbers_mr.csv'   using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/prep_mersenne_primes_mr.csv'    using 1:2 ls 4 title 'Mersenne-Primzahlen'

set output 'plots/prep_frob.tex'
plot 'processed/prep_primes_frob.csv'           using 1:2 ls 1 title 'Primzahlen', \
     'processed/prep_composites_frob.csv'       using 1:2 ls 2 title 'Zusammengesetzte Zahlen', \
     'processed/prep_mersenne_numbers_frob.csv' using 1:2 ls 3 title 'Mersenne-Zahlen', \
     'processed/prep_mersenne_primes_frob.csv'  using 1:2 ls 4 title 'Mersenne-Primzahlen'

# vim: set ft=gnuplot :
