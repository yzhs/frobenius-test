set datafile separator ','

set terminal epslatex color size 15cm,10cm

set xlabel '$\log_2 n$'
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
