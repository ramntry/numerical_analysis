compiler = gcc-4.7

hw09: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c -o $@ -lm

debug : hw09_debug
hw09_debug: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c -o $@ -lm

