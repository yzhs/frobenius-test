# plot_mults.plt -- plot number of multiplications
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
set style line 3 lt rgb "#00B000" lw 2 pt 4
set style line 4 lt rgb "black" lw 2 pt 3
set style line 5 lt rgb "cyan" lw 2
set style line 6 lt rgb "blue" lw 2
set style line 7 lt rgb "violet" lw 2

set xlabel '$\log_2 n$'
set ylabel 'Anzahl der Multiplikationen'

set format x '$2^{%L}$'

set xrange [32:131072]
set logscale x 2
set logscale y

set key top left

C(x) = c*x
D(x) = d*x

fit C(x) 'processed/multiplications_primes.csv' using 1:2 via c
fit D(x) 'processed/multiplications_mersenne_primes.csv' using 1:2 via d

label_C = sprintf("$%d\\log n$", c+0.5)
label_D = sprintf("$%d\\log n$", d+0.5)

set output 'plots/multiplications.tex'
plot 'processed/multiplications_primes.csv'           using 1:2 ls 3 title 'Primzahlen', \
     'processed/multiplications_composites.csv'       using 1:2 ls 4 title 'Zusammengesetzt', \
     'processed/multiplications_mersenne_numbers.csv' using 1:2 ls 5 title 'Mersenne-Zahlen', \
     'processed/multiplications_mersenne_primes.csv'  using 1:2 ls 6 title 'Mersenne-Primzahlen', \
     C(x) title label_C ls 1, \
     D(x) title label_D ls 2

set xrange [2048:131072]
set ylabel 'Laufzeit pro Multiplikation'

# Karatsuba multiplication takes time Θ(x^(log_2 3))
G(x) = g * x ** (log(3)/log(2))

# 3-way Toom-Cook multiplication takes time Θ(x 2^(2√(2log x)) log(x))
H(x) = h * x * 2**(2*sqrt(2*log(x))) * log(x)

fit G(x) 'processed/time_vs_multiplications_primes.csv'          using 3:($2/$1) via g
fit H(x) 'processed/time_vs_multiplications_mersenne_primes.csv' using 3:($2/$1) via h

set output 'plots/time_vs_multiplications.tex'
plot 'processed/time_vs_multiplications_primes.csv'           using 3:($2/$1) ls 3 title 'Primzahlen', \
     'processed/time_vs_multiplications_composites.csv'       using 3:($2/$1) ls 4 title 'Zusammengesetzt', \
     'processed/time_vs_multiplications_mersenne_primes.csv'  using 3:($2/$1) ls 5 title 'Mersenne-Primzahlen', \
     G(x) ls 1 title "Karatsuba", \
     H(x) ls 2 title "Toom-Cook"
     #'processed/time_vs_multiplications_mersenne_numbers.csv' using 3:($2/$1) title 'Mersenne-Zahlen', \

# vim: set ft=gnuplot :
