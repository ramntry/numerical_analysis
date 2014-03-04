compiler = gcc-4.7

hw10: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c -o $@ -lm

debug : hw10_debug
hw10_debug: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c -o $@ -lm

