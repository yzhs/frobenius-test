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

frobenius: frobenius.o helpers.o small_primes.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

frobenius_int: frobenius_int.o helpers_int.o small_primes.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread

miller_rabin: miller_rabin.o helpers.o
	$(CC) $(CFLAGS) -o $@ $^ -lgmp

miller_rabin_int: miller_rabin_int.o helpers_int.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -pthread

run_tests: helpers.o helpers_int.o small_primes.o test/main.o \
	test/test_miller_rabin_long.o test/test_miller_rabin_int.o test/test_frobenius_long.o test/test_frobenius_int.o
	$(CC) $(CFLAGS) -o $@ $^ -lm -lcunit -lgmp

test/%.o: %.c
	$(CC) $(CFLAGS) -DTEST -c -o $@ $<

clean:
	-rm *.o miller_rabin frobenius miller_rabin_int frobenius_int
	-rm test/*.o test/main_int test/main_long run_tests

.PHONY: all clean test test_python
