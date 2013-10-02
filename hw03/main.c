#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "finite_diff_table.h"

double calc(FiniteDiffTable *table, double startX, double step, double point, int isEnd)
{
    size_t const depth = 5;
    assert(depth < table->size);
    size_t offset = (size_t)((point - startX) / step) + isEnd;
    double finiteDiffs[depth];
    for (size_t i = 0; i < depth; ++i) {
        finiteDiffs[i] = fdt_at(*table, offset - i*isEnd, i);
    }
    double const t = (point - (startX + (double)offset * step)) / step;

    printf("\nusing diffs: [ ");
    for (size_t i = 0; i < depth; ++i) {
        printf("%.*lf ", PRINT_NUMERIC_PRECISION, finiteDiffs[i]);
    }
    printf("] and t = %lf\n", t);

    double acc = 0.0;
    for (size_t i = depth - 1; i > 0; --i) {
        acc += finiteDiffs[i];
        double const delta = (double)(i - 1);
        double const new_t = (t - delta + 2.0*delta*isEnd) / (double)i;
        acc *= new_t;
    }
    return acc + finiteDiffs[0];
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input filename>\n", argv[0]);
        return EXIT_FAILURE;
    }
    FILE *inputFile = fopen(argv[1], "r");

    double startX = 0.0;
    double step = 0.0;
    size_t tableSize = 0;
    fscanf(inputFile, "%lf%lf%zu", &startX, &step, &tableSize);

    Vector values;
    initVector(&values, tableSize);
    readVector(inputFile, &values, tableSize);

    double beginX = 0.0;
    double endX = 0.0;
    fscanf(inputFile, "%lf%lf", &beginX, &endX);
    fclose(inputFile);

    FiniteDiffTable table = createFiniteDiffTable(&values);
    disposeVector(&values);

    prettyPrintFiniteDiffTable(&table, stdout);

    double beginValue = calc(&table, startX, step, beginX, 0);
    printf("x = %.*lf\tf(x) = %.*lf\n", PRINT_NUMERIC_PRECISION, beginX, PRINT_NUMERIC_PRECISION, beginValue);

    double endValue = calc(&table, startX, step, endX, 1);
    printf("x = %.*lf\tf(x) = %.*lf\n", PRINT_NUMERIC_PRECISION, endX, PRINT_NUMERIC_PRECISION, endValue);

    disposeFiniteDiffTable(&table);
    return 0;
}

