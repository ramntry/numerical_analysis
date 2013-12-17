compiler = gcc-4.8

hw08: main.c ../support/vector.c ../support/vector.h ../support/polynomial.c ../support/polynomial.h ../hw03/finite_diff_table.c ../hw03/finite_diff_table.h Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c ../support/vector.c ../support/polynomial.c ../hw03/finite_diff_table.c -o $@ -lm

debug : hw08_debug
hw08_debug: main.c ../support/vector.c ../support/vector.h ../support/polynomial.c ../support/polynomial.h ../hw03/finite_diff_table.c ../hw03/finite_diff_table.h Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c ../support/vector.c ../support/polynomial.c ../hw03/finite_diff_table.c -o $@ -lm

