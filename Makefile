CC = clang
ifeq ($(CC),clang)
DEBUG := -DDEBUG -g -Wall -Werror -Wextra -Wmost -Weverything -Wno-pointer-arith -Wno-empty-translation-unit
else
DEBUG := -DDEBUG -g -Wall -Werror -Wextra -Wno-pointer-arith
endif
OPT = -O3 -mtune=native -march=native
CFLAGS = -std=c11 $(DEBUG) $(OPT)


all: miller_rabin miller_rabin_int frobenius frobenius_int

test: run_tests
	./run_tests


frobenius: run_frobenius.o frobenius.o helpers.o small_primes.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

frobenius_int: run_frobenius_int.o frobenius_int.o helpers_int.o small_primes.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread

miller_rabin: run_miller_rabin.o miller_rabin.o helpers.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

miller_rabin_int: run_miller_rabin_int.o miller_rabin_int.o helpers_int.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread


run_tests: helpers.o helpers_int.o small_primes.o test/main.o \
	test/test_miller_rabin_long.o test/test_miller_rabin_int.o test/test_frobenius_long.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -lcunit -lgmp
#test/test_frobenius_int.o


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
	-rm *.o miller_rabin frobenius miller_rabin_int frobenius_int
	-rm test/*.o run_tests

.PHONY: all clean test test_python
