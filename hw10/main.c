#include <complex.h>
#include <stdio.h>
#include <math.h>

#define SIZE 4
#define DOUBLE_WIDTH 8

// 27
complex double matrix[SIZE][SIZE] = {
  { 2.34,  0.91,  0.91, -0.5  },
  { 0.91, -2.12,  0.91,  1.02 },
  { 0.91,  0.91,  2.07,  1.02 },
  {-0.5,   1.02,  1.02,  2.14 }
};

complex double decomposed[SIZE][SIZE];
complex double control[SIZE][SIZE];

complex double vector[SIZE] = {
   5.7565,
  -3.5570,
   3.9305,
   0.3590
};

complex double vector_y[SIZE];
complex double vector_x[SIZE];
complex double control_vector[SIZE];

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

void print_complex(complex double z) {
  printf(" % -7.4f%+.4fi", creal(z), cimag(z));
}

void print_matrix(complex double m[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      print_complex(m[i][j]);
    printf("\n");
  }
  printf("\n");
}

void print_vector(complex double v[SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    print_complex(v[i]);
    printf("\n");
  }
  printf("\n");
}

void back_pass_y()
{
  for (int i = 0; i < SIZE; ++i) {
    vector_y[i] = vector[i];
    for (int j = 0; j < i; ++j)
      vector_y[i] -= decomposed[j][i] * vector_y[j];
    vector_y[i] /= decomposed[i][i];
  }
}

void back_pass_x()
{
  for (int i = SIZE - 1; i >= 0; --i) {
    vector_x[i] = vector_y[i];
    for (int j = i + 1; j < SIZE; ++j)
      vector_x[i] -= decomposed[i][j] * vector_x[j];
    vector_x[i] /= decomposed[i][i];
  }
}

void recalc_vector()
{
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      control_vector[i] += vector_x[j] * matrix[i][j];
  }
}

int main() {
  printf("A =\n");
  print_matrix(matrix);
  printf("b = \n");
  print_vector(vector);

  decompose();
  printf("S =\n");
  print_matrix(decomposed);

  calc_control();
  printf("S'S =\n");
  print_matrix(control);

  back_pass_y();
  printf("y (where S'y = b, y = Sx) = \n");
  print_vector(vector_y);

  back_pass_x();
  printf("x = \n");
  print_vector(vector_x);

  recalc_vector();
  printf("Ax = \n");
  print_vector(control_vector);
  return 0;
}
