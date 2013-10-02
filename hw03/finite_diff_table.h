#ifndef FINITE_DIFF_TABLE_H_
#define FINITE_DIFF_TABLE_H_

#include "../support/vector.h"

typedef struct
{
    Vector triangle;
    size_t size;
}
FiniteDiffTable;

FiniteDiffTable createFiniteDiffTable(Vector const *values);
void disposeFiniteDiffTable(FiniteDiffTable *table);
void printFiniteDiffTable(FiniteDiffTable const *table, FILE *fd);

#define fdt_at(table, i, j) \
    (table).triangle.values[((j)*(2*(table).size - (j) + 1)) / 2 + (i)]

#endif
