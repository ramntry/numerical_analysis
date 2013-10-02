#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "finite_diff_table.h"

double calc(FiniteDiffTable *table, double startX, double step, double point, int isEnd)
{
    size_t const depth = 5;
    assert(depth < table->size);

    size_t const offset = (size_t)((point - startX) / step) + isEnd;
    Vector finiteDiffs;
    initVector(&finiteDiffs, depth);
    for (size_t i = 0; i < depth; ++i) {
        append(&finiteDiffs, fdt_at(*table, offset - i*isEnd, i));
    }
    double const t = (point - (startX + (double)offset * step)) / step;

    printf("\nusing diffs: ");
    printVector(&finiteDiffs, stdout);
    printf(" and t = %lf\n", t);

    double acc = 0.0;
    double sign = 2.0 * (double)isEnd - 1.0;
    for (size_t i = depth - 1; i > 0; --i) {
        acc += finiteDiffs.values[i];
        acc *= (t + sign*(double)(i - 1)) / (double)i;
    }
    acc += finiteDiffs.values[0];

    disposeVector(&finiteDiffs);
    return acc;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input filename>\n", argv[0]);
        return EXIT_FAILURE;
    }
    FILE *inputFile = fopen(argv[1], "r");
    if (inputFile == NULL) {
        fprintf(stderr, "Can not open file `%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    double startX = 0.0;
    double step = 0.0;
    size_t tableSize = 0;
    if (fscanf(inputFile, "%lf%lf%zu", &startX, &step, &tableSize) != 3) {
        fprintf(stderr, "Can not read first x value OR step value OR table size\n");
        return EXIT_FAILURE;
    }

    Vector values;
    initVector(&values, tableSize);
    if (readVector(inputFile, &values, tableSize) != 0) {
        fprintf(stderr, "Can not read %zu floats from table\n", tableSize);
        return EXIT_FAILURE;
    }

    double beginX = 0.0;
    double endX = 0.0;
    if (fscanf(inputFile, "%lf%lf", &beginX, &endX) != 2) {
        fprintf(stderr, "Can not read x values for begin OR end of table\n");
        return EXIT_FAILURE;
    }
    fclose(inputFile);

    FiniteDiffTable table = createFiniteDiffTable(&values);
    disposeVector(&values);

    prettyPrintFiniteDiffTable(&table, stdout);

    double beginValue = calc(&table, startX, step, beginX, /*isEnd=*/0);
    printf("x = %.*lf\tf(x) = %.*lf\n", PRINT_NUMERIC_PRECISION, beginX, PRINT_NUMERIC_PRECISION, beginValue);

    double endValue = calc(&table, startX, step, endX, /*isEnd=*/1);
    printf("x = %.*lf\tf(x) = %.*lf\n", PRINT_NUMERIC_PRECISION, endX, PRINT_NUMERIC_PRECISION, endValue);

    disposeFiniteDiffTable(&table);
    return 0;
}

