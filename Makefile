CC=gcc -std=gnu99
DEBUG=-DDEBUG -g
OPT=

all: miller_rabin frobenius

frobenius: small_primes.c
frobenius_int: small_primes.c

frobenius: frobenius.o helpers.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

miller_rabin: miller_rabin.o helpers.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

frobenius_int: frobenius.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm

miller_rabin_int: miller_rabin.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lm

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
