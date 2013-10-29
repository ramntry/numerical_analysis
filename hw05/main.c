#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "../support/vector.h"

size_t const numof_my_xs = 5;
double pi;
double min_x;
double max_x;

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

double my_f(double x)
{
    return sin(x / 3.0);
}

double my_g(double x)
{
    return x + sin(10.0 * x);
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
        tableAppend(&table, curr_x);
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
        tableAppend(&doubled_table, prev_x + 0.5 * (curr_x - prev_x));
        tableAppend(&doubled_table, curr_x);
        prev_x = curr_x;
    }
    return doubled_table;
}

int main()
{
    pi = 4.0 * atan(1);
    min_x = 0.0;
    max_x = 3.0 * pi;

    Table my_table = createMyTable();
    Table doubled_table = doubleTable(&my_table, min_x);
    Table left_table = createTable(min_x, max_x, 0);
    Table middle_table = createTable(min_x, max_x, 1);
    Table right_table = createTable(min_x, max_x, 2);

    printTable(&my_table, "Initial table");
    printTable(&doubled_table, "Doubled table");
    printTable(&left_table, "Left table");
    printTable(&middle_table, "Middle table");
    printTable(&right_table, "Right table");

    printf("f(%f) = %f\n", pi, lagrangeValue(&my_table, pi));

    disposeTable(&right_table);
    disposeTable(&middle_table);
    disposeTable(&left_table);
    disposeTable(&doubled_table);
    disposeTable(&my_table);
    return 0;
}
