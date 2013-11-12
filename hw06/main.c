#include <stdio.h>
#include <math.h>

double f(double x)
{
  return 1.0 / (0.28 + sinh(x));
}

double const a_lim = 0.0;
double const b_lim = 0.4;
double const max_d2f = 2.0/(0.28*0.28*0.28);
double const max_d4f = (24.9326*(1-0.0789106)*(1-0.0778541)*(1.0778541)*(1.0789106))/(0.28*0.28*0.28*0.28*0.28);

size_t const smallNumSegments = 8;
size_t const bigNumSegments = 16;
size_t const accurateNumSegments = 20*1000;


typedef double (*Function)(double);
typedef double (*Method)(Function, double, double, size_t);
typedef double (*MethodMaxError)(size_t);

double rectangleMethod(Function func, double a, double b, size_t numSegments)
{
  double const step = (b - a) / numSegments;
  double currPoint = a + 0.5 * step;
  double acc = 0.0;
  for (size_t i = 0; i < numSegments; ++i) {
    acc += func(currPoint);
    currPoint += step;
  }
  return step * acc;
}

double rectangleMethodMaxError(size_t numSegments)
{
  double const width = b_lim - a_lim;
  return (width*width*width * max_d2f) / (24.0*numSegments*numSegments);
}

double trapezoidalRule(Function func, double a, double b, size_t numSegments)
{
  double const step = (b - a) / numSegments;
  double currPoint = a + step;
  double acc = 0.5 * (func(a) + func(b));
  for (size_t i = 1; i < numSegments; ++i) {
    acc += func(currPoint);
    currPoint += step;
  }
  return step * acc;
}

double trapezoidalRuleMaxError(size_t numSegments)
{
  return 2.0 * rectangleMethodMaxError(numSegments);
}

double simpsonsRule(Function func, double a, double b, size_t numSegments)
{
  return (trapezoidalRule(func, a, b, numSegments) +
    2.0 * rectangleMethod(func, a, b, numSegments)) / 3.0;
}

double simpsonsRuleMaxError(size_t numSegments)
{
  double const width = b_lim - a_lim;
  double const width2 = width*width;
  double const numSegments2 = numSegments*numSegments;
  return (width2*width2*width * max_d4f) / (2880.0*numSegments2*numSegments2);
}

double accurateByRunge(size_t numSegments)
{
  return (16.0*simpsonsRule(f, a_lim, b_lim, numSegments) -
               simpsonsRule(f, a_lim, b_lim, numSegments/2)) / 15.0;
}

void reportForMethod(Method method, MethodMaxError maxError, double accurate, const char *name)
{
  double const rectSmall = method(f, a_lim, b_lim, smallNumSegments);
  double const rectBig = method(f, a_lim, b_lim, bigNumSegments);
  printf("%-11s | %2zu | %11.9f | %11.9f | %11.9f\n",
      name, smallNumSegments, rectSmall, fabs(accurate - rectSmall), maxError(smallNumSegments));
  printf("            | %2zu | %11.9f | %11.9f | %11.9f\n",
      bigNumSegments, rectBig, fabs(accurate - rectBig), maxError(bigNumSegments));
  printf("------------+----+-------------+-------------+------------\n");
}

int main()
{
  double const accurate = accurateByRunge(accurateNumSegments);
  printf("Method      | N  | S           | Error       | A\n");
  printf("============+====+=============+=============+============\n");
  reportForMethod(rectangleMethod, rectangleMethodMaxError, accurate, "Rectange");
  reportForMethod(trapezoidalRule, trapezoidalRuleMaxError, accurate, "Trapezoidal");
  reportForMethod(simpsonsRule, simpsonsRuleMaxError, accurate, "Simpson's");
  double const runge = accurateByRunge(bigNumSegments);
  printf("%-11s | %2zu | %11.9f | %11.9f | %11.9f\n",
      "Runge", bigNumSegments, runge, fabs(accurate - runge), NAN);
}
