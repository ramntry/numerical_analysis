hw01: main.c ../support/polynomial.h ../support/polynomial.c ../support/vector.h ../support/vector.c Makefile
	cc -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/polynomial.c ../support/vector.c -o $@ -lm

debug : hw01_debug
hw01_debug: main.c ../support/polynomial.h ../support/polynomial.c ../support/vector.h ../support/vector.c Makefile
	cc -std=c99 -pedantic -Wall -g3 main.c ../support/polynomial.c ../support/vector.c -o $@ -lm
