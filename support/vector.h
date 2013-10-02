#ifndef VECTOR_H_
#define VECTOR_H_

#include <stdio.h>

#define PRINT_NUMERIC_PRECISION 6

typedef struct
{
    size_t capacity;
    size_t size;
    double *values;
}
Vector;

void initVector(Vector *vector, size_t capacity);
void disposeVector(Vector *vector);
void append(Vector *vector, double value);
void concat(Vector *dst, Vector const *src);
void readVector(FILE *fd, Vector *vector, size_t num);
void printVector(Vector const *vector, FILE *fd);

#endif
