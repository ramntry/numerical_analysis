#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "../support/vector.h"
#include "../support/polynomial.h"

size_t const numof_my_xs = 5;
size_t const numof_checkpoints = 20;
double pi;
double min_x;
double max_x;
int no_random;

typedef double (*Function)(double x);

typedef struct
{
    Vector xs;
    Vector ys;
    Function f;
}
Table;

void initTable(Table *table, size_t capacity, Function f)
{
    initVector(&table->xs, capacity);
    initVector(&table->ys, capacity);
    table->f = f;
}

void disposeTable(Table *table)
{
    disposeVector(&table->xs);
    disposeVector(&table->ys);
    table->f = NULL;
}

void tableAppend(Table *table, double x)
{
    assert(table->f(x) && "Function stored in the table must be not NULL");
    append(&table->xs, x);
    append(&table->ys, table->f(x));
}

void resetFunction(Table *table, Function new_f)
{
    table->f = new_f;
    for (size_t i = 0; i < table->xs.size; ++i) {
        table->ys.values[i] = new_f(table->xs.values[i]);
    }
}

void printTable(Table const *table, char const *name)
{
    printf("%s:\nxs: ", name);
    printVector(&table->xs, stdout);
    printf("\nys: ");
    printVector(&table->ys, stdout);
    printf("\n\n");
}

double frand()
{
    return no_random ? 0.5 : (double)rand() / RAND_MAX;
}

long fact(long n)
{
    long acc = 1;
    for (long i = n; i > 0; --i) {
        acc *= i;
    }
    return acc;
}

double my_f(double x)
{
    return sin(x / 3.0);
}

double error_coeff_f(long size)
{
    return pow(3.0, -size) / fact(size);
}

double my_g(double x)
{
    return x + sin(10.0 * x);
}

double error_coeff_g(long size)
{
    return pow(10.0, size) / fact(size);
}

double errorUpperbound(Table const *table, double x, double error_coeff)
{
    double acc = 1.0;
    for (size_t j = 0; j < table->xs.size; ++j) {
        acc *= x - table->xs.values[j];
    }
    return fabs(acc * error_coeff);
}

double lagrangeValue(Table const *table, double x)
{
    size_t const size = table->xs.size;
    double acc = 0.0;
    for (size_t k = 0; k < size; ++k) {
        double numerator = 1.0;
        for (size_t j = 0; j < size; ++j) {
            numerator *= (j == k ? 1.0 : x - table->xs.values[j]);
        }
        double denominator = 1.0;
        double const xk = table->xs.values[k];
        for (size_t j = 0; j < size; ++j) {
            denominator *= (j == k ? 1.0 : xk - table->xs.values[j]);
        }
        acc += (numerator * table->ys.values[k]) / denominator;
    }
    return acc;
}

Polynomial createLagrangePolynomial(Table const *table)
{
    size_t const size = table->xs.size;
    Polynomial acc;
    initPolynomial(&acc, size - 1);
    setToScalar(&acc, 0.0);
    Polynomial term;
    initPolynomial(&term, acc.deg);
    for (size_t k = 0; k < size; ++k) {
        setToScalar(&term, 1.0);
        for (size_t j = 0; j < size; ++j) {
            if (j != k) {
                addRoot(&term, table->xs.values[j]);
            }
        }
        double denominator = 1.0;
        double const xk = table->xs.values[k];
        for (size_t j = 0; j < size; ++j) {
            denominator *= (j == k ? 1.0 : xk - table->xs.values[j]);
        }
        multiplyByScalar(&term, table->ys.values[k] / denominator);
        addPolynomial(&acc, &term);
    }
    disposePolynomial(&term);
    return acc;
}

void printReport(Table const *table, char const *name)
{
    putchar('\n');
    printTable(table, name);
    double const error_coeff = (table->f == my_f ? error_coeff_f : error_coeff_g)(table->xs.size);
    char const func_char = table->f == my_f ? 'f' : 'g';
    double const step = (max_x - min_x) / numof_checkpoints;
    double curr_x = min_x;
    double max_error = 0.0;
    double max_upperbound = 0.0;
    Polynomial polynomial = createLagrangePolynomial(table);
    printf("Ln(%c, x) = ", func_char);
    printPolynomial(&polynomial, stdout);
    printf("\n\nNo | x            | %c(x)         | Ln(%c, x)         | error            | A                      | error <= A\n", func_char, func_char);
    printf("---+--------------+--------------+------------------+------------------+------------------------+-----------\n");
    for (size_t i = 1; i <= numof_checkpoints; ++i) {
        double const x = curr_x + frand() * step;
        double const y = table->f(x);
        double const lagrange = lagrangeValue(table, x);

        assert(fabs(calcValue(&polynomial, x) - lagrange) < 1.e-6);

        double const error = fabs(y - lagrange);
        double const error_upperbound = errorUpperbound(table, x, error_coeff);
        printf("%2zu | %+.9f | %+.9f | %+16.9f | %+16.9f | %+22.9f | %s\n"
                , i, x, y, lagrange, error, error_upperbound, error <= error_upperbound ? "Yes" : "No");
        curr_x += step;
        max_error = error > max_error ? error : max_error;
        max_upperbound = error_upperbound > max_upperbound ? error_upperbound : max_upperbound;
    }
    printf("\nmax error = %+.9f\nmax A     = %+.9f\n\n", max_error, max_upperbound);
    disposePolynomial(&polynomial);
}

Table createMyTable()
{
    Table table;
    initTable(&table, numof_my_xs, my_f);
    tableAppend(&table, 0.25 * pi);
    tableAppend(&table, 0.5  * pi);
    tableAppend(&table, 1.   * pi);
    tableAppend(&table, 2.   * pi);
    tableAppend(&table, 3.   * pi);
    return table;
}

Table createTable(double min_x, double max_x, size_t num_part)
{
    size_t const numof_parts = 3;
    double const part_size = (max_x - min_x) / numof_parts;
    double const step = part_size / numof_my_xs;
    Table table;
    initTable(&table, numof_my_xs, my_f);
    double curr_x = min_x + num_part * part_size;
    for (size_t i = 0; i < numof_my_xs; ++i) {
        tableAppend(&table, curr_x + frand() * step);
        curr_x += step;
    }
    return table;
}

Table doubleTable(Table const *table, double min_x)
{
    size_t const old_size = table->xs.size;
    Table doubled_table;
    initTable(&doubled_table, old_size * 2, table->f);
    double prev_x = min_x;
    for (size_t i = 0; i < old_size; ++i) {
        double curr_x = table->xs.values[i];
        assert(curr_x > prev_x && "table must be ascending and starts above min_x");
        tableAppend(&doubled_table, prev_x + (0.25 + 0.5 * frand()) * (curr_x - prev_x));
        tableAppend(&doubled_table, curr_x);
        prev_x = curr_x;
    }
    return doubled_table;
}

int main(int argc, char **argv)
{
#ifndef NDEBUG
    srand(time(NULL));
#endif
    pi = 4.0 * atan(1);
    min_x = 0.0;
    max_x = 3.0 * pi;
    if (argc > 1 && strcmp(argv[1], "-n") == 0) {
        no_random = 1;
    }

    Table my_table = createMyTable();
    Table doubled_table = doubleTable(&my_table, min_x);
    Table left_table = createTable(min_x, max_x, 0);
    Table middle_table = createTable(min_x, max_x, 1);
    Table right_table = createTable(min_x, max_x, 2);

    printReport(&my_table, "Initial table");
    printReport(&doubled_table, "Doubled table");
    printReport(&left_table, "Left table");
    printReport(&middle_table, "Middle table");
    printReport(&right_table, "Right table");

    resetFunction(&right_table, my_g);
    resetFunction(&middle_table, my_g);
    resetFunction(&left_table, my_g);
    resetFunction(&doubled_table, my_g);
    resetFunction(&my_table, my_g);

    printReport(&my_table, "Initial table");
    printReport(&doubled_table, "Doubled table");
    printReport(&left_table, "Left table");
    printReport(&middle_table, "Middle table");
    printReport(&right_table, "Right table");

    disposeTable(&right_table);
    disposeTable(&middle_table);
    disposeTable(&left_table);
    disposeTable(&doubled_table);
    disposeTable(&my_table);
    return 0;
}
