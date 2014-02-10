CC=gcc -std=gnu99
DEBUG=-DDEBUG -g
#OPT=-O3 -mtune=native -march=native -ffast-math -funroll-all-loops
OPT=-O3

all: miller_rabin miller_rabin_int frobenius frobenius_int

frobenius: small_primes.c
frobenius_int: small_primes.c

frobenius: frobenius.o helpers.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

miller_rabin: miller_rabin.o helpers.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

frobenius_int: frobenius_int.o helpers_int.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm -pthread

miller_rabin_int: miller_rabin_int.c helpers_int.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm -pthread

%.o: %.c
	$(CC) $(DEBUG) $(OPT) -c -o $@ $^

.PHONY: all clean test

clean:
	-rm *.o miller_rabin frobenius miller_rabin_int frobenius_int
	cd test && make clean

testpython:
	python -m unittest miller_rabin.py
	python -m unittest frobenius_allinone.py

test:
	cd test && make CC="$(CC)" DEBUG="$(DEBUG)" OPT="$(OPT)" run
