set datafile separator ','

set terminal epslatex color size 17cm,10cm

#set style line 1 lt rgb "red" lw 2
#set style line 2 lt rgb "orange" lw 2
#set style line 3 lt rgb "#f0f000" lw 2
#set style line 4 lt rgb "green" lw 2
#set style line 5 lt rgb "cyan" lw 2
#set style line 6 lt rgb "blue" lw 2
#set style line 7 lt rgb "violet" lw 2

set xlabel '$\lfloor\log_2 n\rfloor$'
set ylabel 'Anzahl der Multiplikationen'

set format x '$2^{%L}$'

set xrange [32:131072]
set logscale x 2
set logscale y

set key top left

C(x) = c*x
D(x) = d*x

fit C(x) 'data/multiplications_primes.csv' using 1:2 via c
fit D(x) 'data/multiplications_mersenne_primes.csv' using 1:2 via d

label_C = sprintf("$%d\\log n$", c+0.5)
label_D = sprintf("$%d\\log n$", d+0.5)

set output 'pic/multiplications.tex'
plot 'data/multiplications_primes.csv'           using 1:2 title 'Primzahlen', \
     'data/multiplications_composites.csv'       using 1:2 title 'Zusammengesetzt', \
     'data/multiplications_mersenne_numbers.csv' using 1:2 title 'Mersenne-Zahlen', \
     'data/multiplications_mersenne_primes.csv'  using 1:2 title 'Mersenne-Primzahlen', \
     C(x) title label_C ls 1, \
     D(x) title label_D ls 2

set xrange [2048:131072]
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
     'data/time_vs_multiplications_mersenne_primes.csv'  using 3:($2/$1) title 'Mersenne-Primzahlen' ls 4, \
     G(x) ls 2 title "Karatsuba", \
     H(x) ls 3 title "Toom-Cook"
     #'data/time_vs_multiplications_mersenne_numbers.csv' using 3:($2/$1) title 'Mersenne-Zahlen', \
