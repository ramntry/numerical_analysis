compiler = gcc-4.8

hw06: main.c ../support/vector.c ../support/vector.h Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/vector.c -o $@ -lm

debug : hw06_debug
hw06_debug: main.c ../support/vector.c ../support/vector.h Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c ../support/vector.c -o $@ -lm
