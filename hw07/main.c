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

typedef double (*Function)(double);

double pi;

double a = 1.e-6;
double b = 0.5;

double alpha = -0.29;

double p(double x)
{
  return pow(x, alpha);
}

double func(double x)
{
  return cos(x);
}

double xs_n2[] = { 0.0790589, 0.381857 };
double xs_n4[] = { 0.0242923, 0.149922, 0.32521, 0.462963 };

double rectangleMethod(Function f, double a, double b, size_t numSegments)
{
  double const step = (b - a) / numSegments;
  double currPoint = a + 0.5 * step;
  double acc = 0.0;
  for (size_t i = 0; i < numSegments; ++i) {
    acc += f(currPoint);
    currPoint += step;
  }
  return step * acc;
}

double trapezoidalRule(Function f, double a, double b, size_t numSegments)
{
  double const step = (b - a) / numSegments;
  double currPoint = a + step;
  double acc = 0.5 * (f(a) + f(b));
  for (size_t i = 1; i < numSegments; ++i) {
    acc += f(currPoint);
    currPoint += step;
  }
  return step * acc;
}

double simpsonsRule(Function f, double a, double b, size_t numSegments)
{
  return (trapezoidalRule(f, a, b, numSegments) +
    2.0 * rectangleMethod(f, a, b, numSegments)) / 3.0;
}

double integrate(Function f)
{
  size_t numSegments = 16;
  double coeff = 16.0;
  return (coeff * simpsonsRule(f, a, b, numSegments) -
                  simpsonsRule(f, a, b, numSegments / 2)) / (coeff - 1.0);
}

double integral_n2_a1_function(double x)
{
  return p(x) * (x - xs_n2[1]) / (xs_n2[0] - xs_n2[1]);
}

double integral_n2_a2_function(double x)
{
  return p(x) * (x - xs_n2[0]) / (xs_n2[1] - xs_n2[0]);
}

double omega_prime(double x)
{
  double a = xs_n4[0];
  double b = xs_n4[1];
  double c = xs_n4[2];
  double d = xs_n4[3];
  return -3*x*x*(a+b+c+d) + 2*x*(a*(b+c+d) + b*(c+d) + c*d) - a*b*c - a*b*d - a*c*d - b*c*d + 4*x*x*x;
}

double integral_n4_a1_function(double x)
{
  double op = omega_prime(xs_n4[0]);
  //printf("%.13f\n", op);
  return p(x) * (x - xs_n4[1]) * (x - xs_n4[2]) * (x - xs_n4[3]) / op;
}

double integral_n4_a2_function(double x)
{
  double op = omega_prime(xs_n4[1]);
  //printf("%.13f\n", op);
  return p(x) * (x - xs_n4[0]) * (x - xs_n4[2]) * (x - xs_n4[3]) / op;
}

double integral_n4_a3_function(double x)
{
  double op = omega_prime(xs_n4[2]);
  //printf("%.13f\n", op);
  return p(x) * (x - xs_n4[0]) * (x - xs_n4[1]) * (x - xs_n4[3]) / op;
}

double integral_n4_a4_function(double x)
{
  double op = omega_prime(xs_n4[3]);
  //printf("%.13f\n", op);
  return p(x) * (x - xs_n4[0]) * (x - xs_n4[1]) * (x - xs_n4[2]) / op;
}

double integral_n2()
{
  double a2 = 0.365517 /*integrate(integral_n2_a1_function)*/;
  double a1 = 0.495498 /*integrate(integral_n2_a2_function)*/;
  printf("\nA = %g, %g", a1, a2);
  return a1 * func(xs_n2[0]) + a2 * func(xs_n2[1]);
}

double integral_n4()
{
  double a1 = 0.22357 /*integrate(integral_n4_a1_function)*/;
  double a2 = 0.284888 /*integrate(integral_n4_a2_function)*/;
  double a3 = 0.236758 /*integrate(integral_n4_a3_function)*/;
  double a4 = 0.115799 /*integrate(integral_n4_a4_function)*/;
  //printf("\nA = %g, %g, %g, %g", a1, a2, a3, a4);
  return a1 * func(xs_n4[0]) + a2 * func(xs_n4[1]) + a3 * func(xs_n4[2]) + a4 * func(xs_n4[3]);
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    pi = 4.0 * atan(1);
    printf("n = 2, x: ");
    for (size_t i = 0; i < 2; ++i)
      printf("%g; ", xs_n2[i]);
    printf("\nS = %g\n", integral_n2());
    printf("\nn = 4, x: ");
    for (size_t i = 0; i < 4; ++i)
      printf("%g; ", xs_n4[i]);
    printf("\nS = %g\n", integral_n4());

    return 0;
}
