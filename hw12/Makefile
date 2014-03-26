compiler = gcc-4.7

hw12: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c -o $@ -lm

debug : hw12_debug
hw12_debug: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c -o $@ -lm

