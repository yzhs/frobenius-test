set datafile separator ','

set terminal epslatex color size 15cm,10cm

set xlabel '$\log_2 n$'
set ylabel 'Anzahl der Multiplikationen'


set xrange [32:65536*2]
set logscale x 2
set logscale y

set key top left

f(x) = (3 + c)*x
g(x) = (3 + d)*x
h(x) = (3 + e)*x

fit f(x) 'data/multiplications_primes.csv' using 1:2 via c
fit g(x) 'data/multiplications_composites.csv' using 1:2 via d
fit h(x) 'data/multiplications_mersenne_primes.csv' using 1:2 via e

label_f = sprintf("$(3 + %d)\\log n$", c+0.5)
label_g = sprintf("$(3 + %d)\\log n$", d+0.5)
label_h = sprintf("$(3 + %d)\\log n$", e+0.5)

set output 'pic/multiplications.tex'
plot 'data/multiplications_primes.csv' using 1:2 title 'Primzahlen', \
     'data/multiplications_composites.csv' using 1:2 title 'Zusammengesetzt', \
     'data/multiplications_mersenne_numbers.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/multiplications_mersenne_primes.csv' using 1:2 title 'Mersenne-Primzahlen', \
     f(x) title label_f ls 1, \
     g(x) title label_g ls 2, \
     h(x) title label_h ls 4, \
     3*x title '$3\log n$' ls 7

set logscale x 10
set xrange [100:1e7]
set xlabel 'Anzahl der Multiplikationen'
set ylabel 'Lauzfeit pro Multiplikation'

set output 'pic/time_vs_multiplications.tex'
plot 'data/time_vs_multiplications_primes.csv' using 1:($2/$1) title 'Primzahlen', \
     'data/time_vs_multiplications_composites.csv' using 1:($2/$1) title 'Zusammengesetzt', \
     'data/time_vs_multiplications_mersenne_numbers.csv' using 1:($2/$1) title 'Mersenne-Zahlen', \
     'data/time_vs_multiplications_mersenne_primes.csv' using 1:($2/$1) title 'Mersenne-Primzahlen', \
