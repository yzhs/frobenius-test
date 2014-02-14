CC=gcc -std=gnu99
DEBUG=-DDEBUG -g -Wall -Werror
#OPT=-O3 -mtune=native -march=native -ffast-math -funroll-all-loops
OPT=-O3

all: miller_rabin miller_rabin_int frobenius frobenius_int

test: run_tests
	./run_tests

frobenius: frobenius.o helpers.o small_primes.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

frobenius_int: frobenius_int.o helpers_int.o small_primes.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm -pthread

miller_rabin: miller_rabin.o helpers.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

miller_rabin_int: miller_rabin_int.o helpers_int.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm -pthread

run_tests: helpers.o helpers_int.o small_primes.o test/main.o \
	test/test_miller_rabin_long.o test/test_miller_rabin_int.o test/test_frobenius_long.o test/test_frobenius_int.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm -lcunit -lgmp

%.o: %.c
	$(CC) $(DEBUG) $(OPT) -c -o $@ $^

test/%.o: %.c
	$(CC) $(DEBUG) $(OPT) -DTEST -c -o $@ $<

clean:
	-rm *.o miller_rabin frobenius miller_rabin_int frobenius_int
	-rm test/*.o test/main_int test/main_long

.PHONY: all clean test test_python
