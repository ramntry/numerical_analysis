#include "finite_diff_table.h"

FiniteDiffTable createFiniteDiffTable(Vector const *values)
{
    FiniteDiffTable table;
    table.size = values->size;
    initVector(&table.triangle, (table.size * (table.size + 1)) / 2);
    concat(&table.triangle, values);
    for (size_t j = 1; j < table.size; ++j) {
        for (size_t i = 0; i < table.size - j; ++i) {
            fdt_at(table, i, j) = fdt_at(table, i + 1, j - 1) - fdt_at(table, i, j - 1);
        }
    }
    return table;
}

void disposeFiniteDiffTable(FiniteDiffTable *table)
{
    disposeVector(&table->triangle);
    table->size = 0;
}

void printFiniteDiffTable(FiniteDiffTable const *table, FILE *fd)
{
    for (size_t i = 0; i < table->size; ++i) {
        for (size_t j = 0; j < table->size - i; ++j) {
            fprintf(fd, "%+.*Lf  ", PRINT_NUMERIC_PRECISION, fdt_at(*table, i, j));
        }
        putchar('\n');
    }
}

void prettyPrintFiniteDiffTable(FiniteDiffTable const *table, FILE *fd)
{
    for (size_t i = 0; i < 2*table->size - 1; ++i) {
        for (size_t j = 0; j < table->size; ++j) {
            if ((i + j) % 2 == 0 && j <= i && i < 2*table->size - j - 1) {
                fprintf(fd, "%+.*Lf  ", PRINT_NUMERIC_PRECISION, fdt_at(*table, (i - j)/2, j));
            } else {
                fprintf(fd, "          ");
            }
        }
        putchar('\n');
    }
}
