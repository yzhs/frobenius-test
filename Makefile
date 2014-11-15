PROFILE=
#CC = clang
#CXX = clang
ifeq ($(CC),clang)
DEBUG = -DDEBUG -g -Wall -Wextra -Wmost -Weverything -Wno-pointer-arith -Wno-empty-translation-unit -Wno-format-nonliteral
else
DEBUG = -DDEBUG -g -Wall -Wextra -Wno-pointer-arith -Wno-format-nonliteral
endif
OPT = -O3 -mtune=native -march=native
CFLAGS = -std=gnu11 $(DEBUG) $(OPT) $(PROFILE)
CXXFLAGS = -std=gnu++11 $(DEBUG) -Wno-c++98-compat-pedantic $(OPT) $(PROFILE)


all: check_all_params check_all_params_long
	#benchmark plots

test: run_tests
	./run_tests

INT_OBJECTS = miller_rabin_int.o frobenius_int.o helpers_int.o
GMP_OBJECTS = miller_rabin.o frobenius.o helpers.o

benchmark: $(INT_OBJECTS) $(GMP_OBJECTS) benchmark.o common.o helpers_int.o helpers.o small_primes.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -lgmp -lm

benchmark.o: benchmark.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $^


SET_TEX = pic/primes.tex pic/composites.tex pic/mersenne_numbers.tex pic/mersenne_primes.tex \
	pic/prep_primes.tex pic/prep_composites.tex pic/prep_mersenne_numbers.tex pic/prep_mersenne_primes.tex
ALG_TEX = pic/gmp.tex pic/mr.tex pic/frob.tex pic/prep_gmp.tex pic/prep_mr.tex pic/prep_frob.tex
MULT_TEX = pic/multiplications.tex
MISC_TEX = pic/all.tex pic/false_positives.tex

plots: $(SET_TEX) $(ALG_TEX) $(MULT_TEX)

data/primes_gmp.csv pic/false_positives.tex: process_timings.jl timings_20140708.csv
	julia $^

$(SET_TEX): plot_sets.plt data/primes_gmp.csv
	gnuplot plot_sets.plt

$(ALG_TEX): plot_algorithms.plt data/primes_gmp.csv
	gnuplot plot_algorithms.plt

$(MULT_TEX): plot_mults.plt data/primes_gmp.csv
	gnuplot plot_mults.plt

frobenius: run_frobenius.o frobenius.o helpers.o small_primes.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

frobenius_int: run_frobenius_int.o frobenius_int.o helpers_int.o small_primes.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread

miller_rabin: run_miller_rabin.o miller_rabin.o helpers.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

miller_rabin_int: run_miller_rabin_int.o miller_rabin_int.o helpers_int.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread

check_all_params_long: check_all_params_long.o frobenius.o helpers.o small_primes.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -lgmp

check_all_params: check_all_params.o frobenius_int.o helpers_int.o small_primes.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lm


run_tests: helpers.o helpers_int.o small_primes.o common.o test/main.o \
	test/test_miller_rabin_long.o test/test_miller_rabin_int.o test/test_frobenius_int.o \
	test/data/primelist.txt
	$(CC) $(CFLAGS) -o $@ $^ -lm -lcunit -lgmp
#	test/test_miller_rabin_long.o test/test_miller_rabin_int.o test/test_frobenius_long.o

# Generate a large list of primes for testing using djb's primegen (http://cr.yp.to/primegen.html)
test/data/primelist.txt:
	primes 2500000000 3000000000 > test/data/primelist.txt

find_non_smooth_numbers: find_non_smooth_numbers.o small_primes.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

nextprime: nextprime.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

add_epsilon_mod_4: add_epsilon_mod_4.o common.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

test/main.o: test/main.c
	$(CC) $(CFLAGS) -c -o $@ $^

test/test_miller_rabin_int.o: test/test_miller_rabin_int.c test/test_miller_rabin_int.h miller_rabin_int.c helpers_int.c helpers_int.h common.h
	$(CC) $(CFLAGS) -c -o $@ test/test_miller_rabin_int.c

test/test_miller_rabin_long.o: test/test_miller_rabin_long.c test/test_miller_rabin_long.h miller_rabin.c helpers.c helpers.h common.h
	$(CC) $(CFLAGS) -c -o $@ test/test_miller_rabin_long.c

test/test_frobenius_int.o: test/test_frobenius_int.c test/test_frobenius_int.h frobenius_int.c helpers_int.c helpers_int.h common.h small_primes.c small_primes.h
	$(CC) $(CFLAGS) -c -o $@ test/test_frobenius_int.c

test/test_frobenius_long.o: test/test_frobenius_long.c test/test_frobenius_long.h frobenius.c helpers.c helpers.h common.h small_primes.c small_primes.h
	$(CC) $(CFLAGS) -c -o $@ test/test_frobenius_long.c


clean:
	-rm *.o test/*.o
	-rm miller_rabin frobenius miller_rabin_int frobenius_int
	-rm run_tests benchmark nextprime find_non_smooth_numbers
	-rm check_all_params check_all_params_long

.PHONY: all clean test test_python
