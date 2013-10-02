#include <assert.h>
#include <stdlib.h>

#include "vector.h"

void initVector(Vector *vector, size_t capacity)
{
    vector->capacity = capacity;
    vector->size = 0;
    if (capacity > 0) {
        vector->values = (double *)malloc(sizeof(double) * capacity);
    } else {
        vector->values = NULL;
    }
}

void disposeVector(Vector *vector)
{
    free(vector->values);
    vector->values = NULL;
    vector->size = 0;
    vector->capacity = 0;
}

void append(Vector *vector, double value)
{
    assert(vector->size < vector->capacity);
    vector->values[vector->size++] = value;
}

void concat(Vector *dst, Vector const *src)
{
    assert(dst->size + src->size < dst->capacity);
    for (size_t i = 0; i < src->size; ++i) {
        dst->values[dst->size + i] = src->values[i];
    }
    dst->size += src->size;
}

void printVector(Vector const *vector, FILE *fd)
{
    fprintf(fd, "[");
    if (vector->size != 0) {
        fprintf(fd, "%.*lf", PRINT_NUMERIC_PRECISION, vector->values[0]);
    }
    for (size_t i = 1; i < vector->size; ++i) {
        fprintf(fd, ", %.*lf", PRINT_NUMERIC_PRECISION, vector->values[i]);
    }
    fprintf(fd, "]");
}

int readVector(FILE *fd, Vector *vector, size_t num)
{
    assert(num <= vector->capacity);
    for (size_t i = 0; i < num; ++i) {
        if (fscanf(fd, "%lf", &vector->values[i]) != 1) {
            return READ_VECTOR_ERROR;
        }
    }
    vector->size = num;
    return 0;
}
