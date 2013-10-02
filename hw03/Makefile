hw03: main.c ../support/vector.c ../support/vector.h finite_diff_table.c finite_diff_table.h Makefile
	cc -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/vector.c finite_diff_table.c -o $@ -lm

debug : hw03_debug
hw03_debug: main.c ../support/vector.c ../support/vector.h finite_diff_table.c finite_diff_table.h  Makefile
	cc -std=c99 -pedantic -Wall -g3 main.c ../support/vector.c finite_diff_table.c -o $@ -lm
