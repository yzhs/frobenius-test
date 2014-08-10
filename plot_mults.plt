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
fit g(x) 'data/multiplications_mersenne_primes.csv' using 1:2 via d
fit h(x) 'data/multiplications_composites.csv' using 1:2 via e

label_f = sprintf("$(3 + %d)\\log n$", c)
label_g = sprintf("$(3 + %d)\\log n$", d)
label_h = sprintf("$(3 + %d)\\log n$", e)

set output 'pic/multiplications.tex'
plot 'data/multiplications_primes.csv' using 1:2 title 'Primzahlen', \
     'data/multiplications_composites.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/multiplications_mersenne_numbers.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/multiplications_mersenne_primes.csv' using 1:2 title 'Mersenne-Primzahlen', \
     f(x) title label_f ls 1, \
     g(x) title label_g ls 4, \
     h(x) title label_h ls 2
