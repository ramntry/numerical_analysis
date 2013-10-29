hw05: main.c ../support/vector.c ../support/vector.h Makefile
	cc -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/vector.c -o $@ -lm

debug : hw05_debug
hw05_debug: main.c ../support/vector.c ../support/vector.h Makefile
	cc -std=c99 -pedantic -Wall -g3 main.c ../support/vector.c -o $@ -lm
