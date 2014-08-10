set datafile separator ','

set terminal epslatex color size 15cm,10cm

set xlabel '$\log_2 n$'
set ylabel 'Anzahl der Multiplikationen'


set xrange [32:100000]
set logscale x 2
set logscale y

set key top left
set ytics 1e-6,1e2,1e8

set output 'pic/multiplications.tex'
plot 'data/multiplications_primes.csv' using 1:2 title 'Primzahlen', \
     'data/multiplications_composites.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/multiplications_mersenne_numbers.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/multiplications_mersenne_primes.csv' using 1:2 title 'Mersenne-Primzahlen'
