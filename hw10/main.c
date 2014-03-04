#include <complex.h>
#include <stdio.h>
#include <math.h>

#define SIZE 5
#define DOUBLE_WIDTH 8

double const eps = 1.0e-5;

// 27
complex double matrix[SIZE][SIZE] = {
  {  1,  3, -2,  0, -2 },
  {  3,  4, -5,  1, -3 },
  { -2, -5,  3, -2,  2 },
  {  0,  1, -2,  5,  3 },
  { -2, -3,  2,  3,  4 }
};

complex double decomposed[SIZE][SIZE];
complex double control[SIZE][SIZE];

complex double vector[SIZE] = {
  0.5,
  5.4,
  5.0,
  7.5,
  3.3
};

void decompose() {
  for (int i = 0; i < SIZE; ++i)
    for (int j = i; j < SIZE; ++j) {
      decomposed[i][j] = matrix[i][j];
      for (int k = 0; k < i; ++k)
        decomposed[i][j] -= decomposed[k][i] * decomposed[k][j];
      if (i == j)
        decomposed[i][j] = csqrt(decomposed[i][j]);
      else
        decomposed[i][j] /= decomposed[i][i];
    }
}

void calc_control() {
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      for (int k = 0; k < SIZE; ++k)
        control[i][j] += decomposed[k][i] * decomposed[k][j];
}

void print_matrix(complex double m[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      printf(" % -4.1f%+.2fi", creal(m[i][j]), cimag(m[i][j]));
    printf("\n");
  }
  printf("\n");
}

int main() {
  print_matrix(matrix);
  decompose();
  print_matrix(decomposed);
  calc_control();
  print_matrix(control);
  return 0;
}
