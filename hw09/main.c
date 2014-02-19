#include <stdio.h>
#include <math.h>

#define SIZE 5
#define DOUBLE_WIDTH 8

// IX

double const eps = 1.0e-5;

double matrix[SIZE][SIZE] = {
  {  5.11, -2.32,  0.46,  1.52,  0.00 },
  {  1.17, -4.08, -3.25,  3.25,  2.34 },
  {  0.00,  1.74,  3.78,  2.15,  1.83 },
  {  0.00,  2.34,  2.05, -3.50,  1.25 },
  {  0.00,  0.00,  1.82,  2.67, -4.79 }
};

double vector[SIZE] = {
  -10.172,
    1.356,
   -1.015,
   -2.295,
    1.425
};

double orig_matrix[SIZE][SIZE];
double control[SIZE];

void calc_control()
{
  for (int i = 0; i < SIZE; ++i) {
    control[i] = vector[i];
    for (int j = 0; j < SIZE; ++j)
      control[i] += matrix[i][j];
  }
}

void print_all(char const *matrix_name, char const *vector1, char const *vector2)
{
  printf(" %-*s | %-*s | %-*s\n", DOUBLE_WIDTH * SIZE, matrix_name, DOUBLE_WIDTH, vector1, DOUBLE_WIDTH, vector2);
  for (int i = 0; i < DOUBLE_WIDTH * (SIZE + 3); ++i)
    putchar('-');
  putchar('\n');
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      printf("% *.3f", DOUBLE_WIDTH, matrix[i][j]);
    printf("  | % *.3f | % *.3f\n", DOUBLE_WIDTH, vector[i], DOUBLE_WIDTH, control[i]);
  }
  putchar('\n');
}

void swap(double *a, double *b)
{
  double const tmp = *a;
  *a = *b;
  *b = tmp;
}

void gaussian()
{
  for (int i = 0; i < SIZE; ++i)
    for (int j = i; j < SIZE; ++j)
      if (fabs(matrix[j][i]) > eps) {
        for (int k = i; k < SIZE; ++k)
          swap(&matrix[i][k], &matrix[j][k]);
        swap(&vector[i], &vector[j]);
        swap(&control[i], &control[j]);

        for (int k = i + 1; k < SIZE; ++k)
          matrix[i][k] /= matrix[i][i];
        vector[i] /= matrix[i][i];
        control[i] /= matrix[i][i];
        matrix[i][i] = 1.0;

        for (int k = i + 1; k < SIZE; ++k) {
          for (int l = i + 1; l < SIZE; ++l)
            matrix[k][l] -= matrix[k][i] * matrix[i][l];
          vector[k] -= matrix[k][i] * vector[i];
          control[k] -= matrix[k][i] * control[i];
          matrix[k][i] = 0.0;
        }
        break;
      }
}

void back_pass()
{
  for (int i = SIZE - 1; i >= 0; --i) {
    control[i] = vector[i];
    for (int j = i + 1; j < SIZE; ++j)
      control[i] -= matrix[i][j] * control[j];
  }
}

void save_orig_matrix()
{
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      orig_matrix[i][j] = matrix[i][j];
}

void recalc_vector()
{
  for (int i = 0; i < SIZE; ++i) {
    vector[i] = 0.0;
    for (int j = 0; j < SIZE; ++j)
      vector[i] += control[j] * orig_matrix[i][j];
  }
}

int main()
{
  save_orig_matrix();
  calc_control();
  print_all("A", "b", "b^");

  gaussian();
  print_all("A~", "b~", "b^~ (exp)");
  calc_control();
  print_all("A~", "b~", "b~^ (act)");

  back_pass();
  print_all("A~", "b~", "x");
  recalc_vector();
  print_all("A~", "b = Ax", "x");
}
