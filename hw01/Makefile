hw01: main.c
	cc -std=c99 -pedantic -Wall -O3 -DNDEBUG $^ -o $@ -lm

debug : hw01_debug
hw01_debug: main.c
	cc -std=c99 -pedantic -Wall -g3 $^ -o $@ -lm
