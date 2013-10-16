hw04: main.c ../support/vector.c ../support/vector.h ../hw03/finite_diff_table.c ../hw03/finite_diff_table.h Makefile
	cc -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/vector.c ../hw03/finite_diff_table.c -o $@ -lm

debug : hw04_debug
hw04_debug: main.c ../support/vector.c ../support/vector.h ../hw03/finite_diff_table.c ../hw03/finite_diff_table.h  Makefile
	cc -std=c99 -pedantic -Wall -g3 main.c ../support/vector.c ../hw03/finite_diff_table.c -o $@ -lm
