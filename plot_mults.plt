set datafile separator ','

set terminal epslatex color size 17cm,10cm

set xlabel '$\lfloor\log_2 n\rfloor$'
set ylabel 'Anzahl der Multiplikationen'


set xrange [32:65536*2]
set logscale x 2
set logscale y

set key top left

C(x) = (3 + c)*x
D(x) = (3 + d)*x
E(x) = (3 + e)*x

fit C(x) 'data/multiplications_primes.csv' using 1:2 via c
fit D(x) 'data/multiplications_composites.csv' using 1:2 via d
fit E(x) 'data/multiplications_mersenne_primes.csv' using 1:2 via e

label_C = sprintf("$(3 + %d)\\log n$", c+0.5)
label_D = sprintf("$(3 + %d)\\log n$", d+0.5)
label_E = sprintf("$(3 + %d)\\log n$", e+0.5)

set output 'pic/multiplications.tex'
plot 'data/multiplications_primes.csv' using 1:2 title 'Primzahlen', \
     'data/multiplications_composites.csv' using 1:2 title 'Zusammengesetzt', \
     'data/multiplications_mersenne_numbers.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/multiplications_mersenne_primes.csv' using 1:2 title 'Mersenne-Primzahlen', \
     C(x) title label_C ls 1, \
     D(x) title label_D ls 2, \
     E(x) title label_E ls 4, \
     3*x title '$3\log n$' ls 7

set xrange [2048:65536*2]
set ylabel 'Laufzeit pro Multiplikation'

# Karatsuba multiplication takes time Θ(x^(log_2 3))
G(x) = g * x ** (log(3)/log(2))

# 3-way Toom-Cook multiplication takes time Θ(x 2^(2√(2log x)) log(x))
H(x) = h * x * 2**(2*sqrt(2*log(x))) * log(x)

fit G(x) 'data/time_vs_multiplications_primes.csv'          using 3:($2/$1) via g
fit H(x) 'data/time_vs_multiplications_mersenne_primes.csv' using 3:($2/$1) via h

set output 'pic/time_vs_multiplications.tex'
plot 'data/time_vs_multiplications_primes.csv'           using 3:($2/$1) title 'Primzahlen', \
     'data/time_vs_multiplications_composites.csv'       using 3:($2/$1) title 'Zusammengesetzt', \
     'data/time_vs_multiplications_mersenne_numbers.csv' using 3:($2/$1) title 'Mersenne-Zahlen', \
     'data/time_vs_multiplications_mersenne_primes.csv'  using 3:($2/$1) title 'Mersenne-Primzahlen', \
     G(x) ls 2 title "Karatsuba", \
     H(x) ls 3 title "Toom-Cook"
