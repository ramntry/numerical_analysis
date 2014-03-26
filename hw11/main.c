#include <complex.h>
#include <stdio.h>
#include <math.h>

#define SIZE 3
#define DOUBLE_WIDTH 8
#define EPS 1.e-8
#define DROP_VALUE 0.1
#define STEP 0.05

void print_vector(double x[SIZE]) {
  for (int i = 0; i < SIZE; ++i)
    printf(" % .*f\n", DOUBLE_WIDTH, x[i]);
}

void print_matrix(double m[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      printf(" % .*f", DOUBLE_WIDTH, m[i][j]);
    printf("\n");
  }
}

// II
double origin[SIZE][SIZE] = {
  {  7.352720,  0.8825500, -2.2700524  },
  {  0.882550,  5.5835140,  0.52816684 },
  { -2.2700524, 0.52816684, 4.4303291  },
};

double matrix[SIZE][SIZE];

double vector[SIZE] = {
  -0.12225372,
  -0.8699978,
   0.9964601
};

double b[SIZE];

int iteration(double m[SIZE][SIZE], double x[SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      if (i != j)
        m[i][j] /= m[i][i];
    b[i] = vector[i] / m[i][i];
    m[i][i] = 0.0;
    x[i] = b[i];
  }
  int iterate = 1;
  int counter = 0;
  while (iterate) {
    double next_x[SIZE];
    iterate = 0;
    for (int i = 0; i < SIZE; ++i) {
      next_x[i] = b[i];
      for (int j = 0; j < SIZE; ++j)
        next_x[i] -= m[i][j] * x[j];
      iterate = fabs(x[i] - next_x[i]) > EPS;
    }
    for (int i = 0; i < SIZE; ++i)
      x[i] = next_x[i];
    ++counter;
  }
  return counter;
}

int nekrasov(double m[SIZE][SIZE], double x[SIZE]) {
  for (int i = 0; i < SIZE; ++i)
    x[i] = b[i];
  int iterate = 1;
  int counter = 0;
  while (iterate) {
    iterate = 0;
    for (int i = 0; i < SIZE; ++i) {
      double next_xi = b[i];
      for (int j = 0; j < SIZE; ++j)
        next_xi -= m[i][j] * x[j];
      iterate = fabs(x[i] - next_xi) > EPS;
      x[i] = next_xi;
    }
    ++counter;
  }
  return counter;
}

void calc_error(double m[SIZE][SIZE], double x[SIZE], double err[SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    err[i] = 0.0;
    for (int j = 0; j < SIZE; ++j)
      err[i] += m[i][j] * x[j];
    err[i] = fabs(err[i] - vector[i]);
  }
}

void test() {
  printf("A =\n");
  print_matrix(origin);
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      matrix[i][j] = origin[i][j];
  double x[SIZE];
  printf("\n== Iteration ==\n");
  int counter = iteration(matrix, x);
  printf("x =\n");
  print_vector(x);
  printf("Number of iterations: %d\n|Ax - b| =\n", counter);
  double err[SIZE];
  calc_error(origin, x, err);
  print_vector(err);

  printf("\n== Nekrasov ==\n");
  counter = nekrasov(matrix, x);
  printf("x =\n");
  print_vector(x);
  printf("Number of iterations: %d\n|Ax - b| =\n", counter);
  calc_error(origin, x, err);
  print_vector(err);
  printf("---------------------------------------\n\n");
}

void raise_force() {
  for (int i = 0; i < SIZE; ++i)
    origin[i][i] *= 10;
}

void back_force() {
  for (int i = 0; i < SIZE; ++i)
    origin[i][i] /= 100;
}

void drop_force(int i, double step) {
  double sum = 0.0;
  for (int j = 0; j < SIZE; ++j)
    if (i != j)
      sum += fabs(origin[i][j]);
  origin[i][i] = sum - step * DROP_VALUE;
}

double find_third() {
  int s = 0;
  for (;; ++s) {
    drop_force(2, -s * STEP);
    for (int i = 0; i < SIZE; ++i)
      for (int j = 0; j < SIZE; ++j)
        matrix[i][j] = origin[i][j];
    double x[SIZE];
    iteration(matrix, x);
    for (int i = 0; i < SIZE; ++i)
      if (isinf(x[i]) || isnan(x[i]))
        break;
      else
        return s * STEP * DROP_VALUE;
  }
  return s * STEP * DROP_VALUE;
}

int main() {
  test();
  raise_force();
  test();
  raise_force();
  test();
  back_force();
  drop_force(0, 1);
  test();
  drop_force(1, 1);
  test();
  double f = find_third();
  test();
  printf("Force = %f\n", f);
  return 0;
}
