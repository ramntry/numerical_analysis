compiler = gcc-4.7

hw11: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -O3 -DNDEBUG main.c -o $@ -lm

debug : hw11_debug
hw11_debug: main.c Makefile
	$(compiler) -std=c99 -pedantic -Wall -g3 main.c -o $@ -lm

