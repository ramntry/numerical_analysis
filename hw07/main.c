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
  size_t numof_xs = (size_t)((length + eps) / h_value) + 1;
  printf("Number of points = %zu\n", numof_xs);
  Vector ys;
  initVector(&ys, numof_xs);
  ys.values[0] = y0;
  for (size_t i = 1; i < numof_xs; ++i) {
    ys.values[i] = ys.values[i - 1] + h_value * f(min_x + i * h_value, ys.values[i - 1]);
  }
  ys.size = numof_xs;
  return ys;
}

void print_euler(double h_value)
{
  printf("Euler with h = %g\n", h_value);
  Vector euler_ys = euler(h_value);
  printVector(&euler_ys, stdout);
  putchar('\n');
  disposeVector(&euler_ys);
}

void print_results()
{
  print_euler(2.0 * h);
  print_euler(1.0 * h);
  print_euler(0.5 * h);
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    pi = 4.0 * atan(1);
    length = max_x - min_x;

    print_results();

    return 0;
}
