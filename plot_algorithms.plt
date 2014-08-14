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


set xrange [32:65536*2]
set logscale x 2
set logscale y

set key top left
set ytics 1e-6,1e2,1e8

# Full test
set output 'pic/gmp.tex'
plot 'data/primes_gmp.csv' using 1:2 title 'Primzahlen', \
     'data/composites_gmp.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/mersenne_numbers_gmp.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/mersenne_primes_gmp.csv' using 1:2 title 'Mersenne-Primzahlen'

set output 'pic/mr.tex'
plot 'data/primes_mr.csv' using 1:2 title 'Primzahlen', \
     'data/composites_mr.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/mersenne_numbers_mr.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/mersenne_primes_mr.csv' using 1:2 title 'Mersenne-Primzahlen'

set output 'pic/frob.tex'
plot 'data/primes_frob.csv' using 1:2 title 'Primzahlen', \
     'data/composites_frob.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/mersenne_numbers_frob.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/mersenne_primes_frob.csv' using 1:2 title 'Mersenne-Primzahlen'


# Just the precomputation
set output 'pic/prep_gmp.tex'
plot 'data/prep_primes_gmp.csv' using 1:2 title 'Primzahlen', \
     'data/prep_composites_gmp.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/prep_mersenne_numbers_gmp.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/prep_mersenne_primes_gmp.csv' using 1:2 title 'Mersenne-Primzahlen'

set output 'pic/prep_mr.tex'
plot 'data/prep_primes_mr.csv' using 1:2 title 'Primzahlen', \
     'data/prep_composites_mr.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/prep_mersenne_numbers_mr.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/prep_mersenne_primes_mr.csv' using 1:2 title 'Mersenne-Primzahlen'

set output 'pic/prep_frob.tex'
plot 'data/prep_primes_frob.csv' using 1:2 title 'Primzahlen', \
     'data/prep_composites_frob.csv' using 1:2 title 'Zusammengesetzte Zahlen', \
     'data/prep_mersenne_numbers_frob.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/prep_mersenne_primes_frob.csv' using 1:2 title 'Mersenne-Primzahlen'
