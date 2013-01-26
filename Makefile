CC=gcc -std=gnu99
DEBUG=-DDEBUG -g
OPT=

all: miller_rabin

miller_rabin: miller_rabin.o
	$(CC) $(DEBUG) $(OPT) -o $@ $^ -lgmp

%.o: %.c
	$(CC) $(DEBUG) $(OPT) -c -o $@ $^

.PHONY: all clean

clean:
	rm miller_rabin miller_rabin.o

test:
	python -m unittest polynomial.py
	python -m unittest miller_rabin.py
	python -m unittest frobenius.py
