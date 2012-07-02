CC=gcc -std=gnu99
DEBUG=-DDEBUG -g

all: miller_rabin.o

miller_rabin: miller_rabin.o
	$(CC) $(DEBUG) -o $@ $^ -lgmp

%.o: %.c
	$(CC) $(DEBUG) -c -o $@ $^

.PHONY: all clean

clean:
	rm miller_rabin.o
