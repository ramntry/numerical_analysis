#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../hw03/finite_diff_table.h"

#define EPS (1.e-6)

double uniCalc(FiniteDiffTable *table, int pos, int offset, double t)
{
    int const depth = 5;
    assert(depth < table->size);

    Vector finiteDiffs;
    initVector(&finiteDiffs, depth);
    for (int i = 0; i < depth; ++i) {
        int secondOffset = offset + ((pos < 2) ? -i*pos : -i/2);
        append(&finiteDiffs, (i < table->size - secondOffset) ? fdt_at(*table, secondOffset, i) : 0.0);
    }

    printf("[%s of table] using diffs: ", pos ? (pos == 1 ? "end" : "middle") : "begin");
    printVector(&finiteDiffs, stdout);
    printf(" and t = %lf\n", t);

    double acc = 0.0;
    for (size_t i = depth - 1; i > 0; --i) {
        int offset = (pos < 2) ? (2*pos - 1)*(i - 1) : i/2 * (i % 2 ? 1 : -1);
        acc += finiteDiffs.values[i];
        acc *= (t + offset) / (double)i;
    }
    acc += finiteDiffs.values[0];

    disposeVector(&finiteDiffs);
    return acc;
}

double calc(FiniteDiffTable *table, double startX, double step, double point, int pos /*0 - begin, 1 - end, 2 - middle*/)
{
    int const offset = (point - startX) / step + pos % 2;
    double const t = (point - (startX + (double)offset * step)) / step;
    return uniCalc(table, pos, offset, t);
}

double phi(FiniteDiffTable *table, double pointY, int pos, int offset, double t)
{
    double const polyValue = uniCalc(table, pos, offset, t);
    double const linKoeff = fdt_at(*table, offset - pos % 2, 1);
    return (pointY - polyValue) / linKoeff + t;
}

double iterate(FiniteDiffTable *table, double startX, double step, double pointY, int pos, int offset)
{
    double prevT = 0.0;
    printf("1: ");
    double currT = phi(table, pointY, pos, offset, prevT);
    for (int i = 2; fabs(prevT - currT) > EPS; ++i) {
        printf("%d: ", i);
        prevT = currT;
        currT = phi(table, pointY, pos, offset, prevT);
    }
    return startX + offset * step + currT * step;
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

    double middleX = 0.0;
    if (fscanf(inputFile, "%lf", &middleX) != 1) {
        fprintf(stderr, "Can not read x values for middle of table\n");
        return EXIT_FAILURE;
    }
    fclose(inputFile);

    FiniteDiffTable table = createFiniteDiffTable(&values);
    disposeVector(&values);

    prettyPrintFiniteDiffTable(&table, stdout);
    putchar('\n');

    double beginValue = calc(&table, startX, step, middleX, /*pos=begin*/0);
    printf("x = %.*lf\tf(x) = %.*lf\n\n", PRINT_NUMERIC_PRECISION, middleX, PRINT_NUMERIC_PRECISION + 1, beginValue);

    double endValue = calc(&table, startX, step, middleX, /*pos=end*/1);
    printf("x = %.*lf\tf(x) = %.*lf\n\n", PRINT_NUMERIC_PRECISION, middleX, PRINT_NUMERIC_PRECISION + 1, endValue);

    double middleValue = calc(&table, startX, step, middleX, /*pos=end*/2);
    printf("x = %.*lf\tf(x) = %.*lf\n\n", PRINT_NUMERIC_PRECISION, middleX, PRINT_NUMERIC_PRECISION + 1, middleValue);

    int const reverseOffset = 4;
    int const reversePos = 2;
    double const reverseY = 1.782125;
    double const reverseResult = iterate(&table, startX, step, reverseY, reversePos, reverseOffset);
    printf("\n[reverse] f(x) = %.*lf\tx = %.*lf\n", PRINT_NUMERIC_PRECISION, reverseY, PRINT_NUMERIC_PRECISION, reverseResult);

    disposeFiniteDiffTable(&table);
    return 0;
}

