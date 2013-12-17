#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "../support/vector.h"
#include "../support/polynomial.h"
#include "../hw03/finite_diff_table.h"

double pi;
double eps = 1.e-6;

double h = 0.1;
double min_x = 0.0;
double max_x = 1.0;
double length = 0.0;

size_t get_numof_xs(double h_value)
{
  return (size_t)((length + eps) / h_value) + 1;
}

double get_x(size_t k, double h_value)
{
  return min_x + k * h_value;
}

double min_a = 14.0;
double max_a = 38.0;
double step_a = 1.0;

double alpha(double a)
{
  return 2.5 + a / 40.0;
}

double current_alpha = 0.0;

double f(double x, double y)
{
  return -y*y + current_alpha * x / (1.0 + x*x);
}

double y0 = -0.4122;

Vector euler(double h_value)
{
  size_t numof_xs = get_numof_xs(h_value);
  Vector ys;
  initVector(&ys, numof_xs);
  ys.values[0] = y0;
  for (size_t i = 1; i < numof_xs; ++i) {
    ys.values[i] = ys.values[i - 1] + h_value * f(get_x(i - 1, h_value), ys.values[i - 1]);
  }
  ys.size = numof_xs;
  return ys;
}

void print_euler(double h_value)
{
  printf("Euler with h = %g\n", h_value);
  Vector euler_ys = euler(h_value);
  printVector(&euler_ys, stdout);
  printf("\n\n");
  disposeVector(&euler_ys);
}

double runge_next_yk(double xk, double yk)
{
  double k1 = h * f(xk, yk);
  double k2 = h * f(xk + 0.5 * h, yk + 0.5 * k1);
  double k3 = h * f(xk + 0.5 * h, yk + 0.5 * k2);
  double k4 = h * f(xk + h, yk + k3);
  return yk + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

Vector runge()
{
  size_t numof_xs = get_numof_xs(h);
  Vector ys;
  initVector(&ys, numof_xs);
  ys.values[0] = y0;
  for (size_t i = 1; i < numof_xs; ++i) {
    ys.values[i] = runge_next_yk(get_x(i - 1, h), ys.values[i - 1]);
  }
  ys.size = numof_xs;
  return ys;
}

void print_runge()
{
  printf("Runge-Kutta\n");
  Vector runge_ys = runge();
  printVector(&runge_ys, stdout);
  printf("\n\n");
  disposeVector(&runge_ys);
}

void printTable(FiniteDiffTable const *table, Vector const *ys, FILE *fd)
{
  size_t row_i = 0;
  int raw_begin = 1;
  for (size_t i = 0; i < 2*table->size - 1; ++i) {
    for (size_t j = 0; j < table->size; ++j) {
      if ((i + j) % 2 == 0 && j <= i && i < 2*table->size - j - 1) {
        if (raw_begin) {
          double x = get_x(row_i, h);
          double y = ys->values[row_i];
          fprintf(fd, "%+.*f  %+.*f  ", PRINT_NUMERIC_PRECISION, x, PRINT_NUMERIC_PRECISION, y);
          ++row_i;
          raw_begin = 0;
        }
        fprintf(fd, "%+.*Lf  ", PRINT_NUMERIC_PRECISION, fdt_at(*table, (i - j)/2, j));
      } else {
        if (raw_begin) {
          fprintf(fd, "                                 ");
          raw_begin = 0;
        } else
          fprintf(fd, "          ");
      }
    }
    putchar('\n');
    raw_begin = 1;
  }
}

double eta(size_t k, double yk)
{
  return h * f(get_x(k, h), yk);
}

size_t const adams_acc = 5;

double adams_next_yk(double yk, FiniteDiffTable const *table)
{
  return fdt_at(*table, table->size - 1 - 0, 0)
       + fdt_at(*table, table->size - 1 - 1, 1) * (  1.0 /   2.0)
       + fdt_at(*table, table->size - 1 - 2, 2) * (  5.0 /  12.0)
       + fdt_at(*table, table->size - 1 - 3, 3) * (  3.0 /   8.0)
       + fdt_at(*table, table->size - 1 - 4, 4) * (251.0 / 720.0)
       + yk;
}

Vector adams()
{
  size_t numof_xs = get_numof_xs(h);
  Vector runge_ys = runge();
  Vector adams_ys;
  initVector(&adams_ys, numof_xs);
  for (size_t i = 0; i < adams_acc; ++i)
    append(&adams_ys, runge_ys.values[i]);
  disposeVector(&runge_ys);

  Vector etas;
  initVector(&etas, numof_xs);
  for (size_t i = 0; i < adams_acc; ++i)
    append(&etas, eta(i, adams_ys.values[i]));
  FiniteDiffTable table = createFiniteDiffTable(&etas);

  for (size_t k = 5; k < numof_xs; ++k) {
    double yk = adams_next_yk(adams_ys.values[k - 1], &table);
    append(&adams_ys, yk);
    append(&etas, eta(k, yk));
    disposeFiniteDiffTable(&table);
    table = createFiniteDiffTable(&etas);
  }

  printTable(&table, &adams_ys, stdout);
  disposeVector(&etas);
  disposeFiniteDiffTable(&table);
  return adams_ys;
}

void print_adams()
{
  printf("Adams\n");
  printf(" x          y          eta       \n");
  Vector adams_ys = adams();
  printVector(&adams_ys, stdout);
  printf("\n\n");
  disposeVector(&adams_ys);
}

void print_results()
{
  print_euler(2.0 * h);
  print_euler(1.0 * h);
  print_euler(0.5 * h);
  print_runge();
  print_adams();
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    pi = 4.0 * atan(1);
    length = max_x - min_x;

    print_results();

    return 0;
}
