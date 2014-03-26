#include <assert.h>
#include <complex.h>
#include <stdio.h>
#include <math.h>

#define SIZE 3
#define DOUBLE_PREC 6
#define EPS 1.e-5

void print_vector(double x[SIZE]) {
  for (int i = 0; i < SIZE; ++i)
    printf(" % .*f\n", DOUBLE_PREC, x[i]);
}

void print_matrix(double m[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j)
      printf(" % .*f", DOUBLE_PREC, m[i][j]);
    printf("\n");
  }
}

// 2
double matrix[SIZE][SIZE] = {
  { 2.23, 1.46, 0.94 },
  { 1.46, 2.75, 1.24 },
  { 0.94, 1.24, 3.67 }
};

double vector[SIZE] = {
  -2.180,
  -0.570,
  -5.205
};

double norm(double x[SIZE]) {
  double acc = 0;
  for (int i = 0; i < SIZE; ++i)
    acc += x[i] * x[i];
  return sqrt(acc);
}

void normalize_by_first(double x[][SIZE], int num) {
  double const norm_coeff = 1 / norm(x[0]);
  for (int k = 0; k < num; ++k)
    for (int i = 0; i < SIZE; ++i)
      x[k][i] *= norm_coeff;
}

double iteration(int *iter_counter) {
  *iter_counter = -1;
  double ys[3][SIZE] = { { 0 } };
  for (int i = 0; i < SIZE; ++i)
    ys[0][i] = ys[1][i] = vector[i];
  for (int k = 1;; ++k) {
    if (k % 3 == 0)
      normalize_by_first(ys, 3);
    double *const curr = ys[ k      % 3];
    double *const next = ys[(k + 1) % 3];
    double *const prev = ys[(k - 1) % 3];
    for (int i = 0; i < SIZE; ++i) {
      next[i] = 0;
      for (int j = 0; j < SIZE; ++j)
        next[i] += matrix[i][j] * curr[j];
    }
    printf("\n%d:\n", k);
    print_vector(curr);
    printf("eigenvalue = %.*f\n", DOUBLE_PREC, next[0] / curr[0]);
    int done = 1;
    for (int i = 0; i < SIZE; ++i)
      if (fabs(next[i] / curr[i] - curr[i] / prev[i]) > EPS) {
        done = 0;
        break;
      }
    if (done) {
      printf("\nIteration: eigenvector, corresponding to max by absolute value eigenvalue:\n");
      print_vector(curr);
      *iter_counter = k;
      return next[0] / curr[0];
    }
  }
}

void jacobi_rotation(double m[SIZE][SIZE]) {
  double max_nondiag = 0.0;
  int max_i = 0;
  int max_j = 0;
  for (int i = 0; i < SIZE; ++i)
    for (int j = i + 1; j < SIZE; ++j)
      if (fabs(m[i][j]) > max_nondiag) {
        max_i = i;
        max_j = j;
        max_nondiag = fabs(m[i][j]);
      }
  double const cur_ii = m[max_i][max_i];
  double const cur_jj = m[max_j][max_j];
  double const theta = 0.5 * atan(2 * m[max_i][max_j] / (cur_ii - cur_jj));
  double const c = cos(theta);
  double const s = sin(theta);
  double row_i[SIZE];
  double row_j[SIZE];
  for (int k = 0; k < SIZE; ++k) {
    row_i[k] = c * m[k][max_i] + s * m[k][max_j];
    row_j[k] = c * m[k][max_j] - s * m[k][max_i];
  }
  double const new_ii = c*c*cur_ii + 2*c*s*max_nondiag + s*s*cur_jj;
  double const new_jj = s*s*cur_ii - 2*c*s*max_nondiag + c*c*cur_jj;
  assert(fabs((c*c - s*s)*max_nondiag - c*s*(cur_ii - cur_jj)) < EPS);
  for (int k = 0; k < SIZE; ++k) {
    m[k][max_i] = m[max_i][k] = row_i[k];
    m[k][max_j] = m[max_j][k] = row_j[k];
  }
  m[max_i][max_i] = new_ii;
  m[max_j][max_j] = new_jj;
  m[max_i][max_j] = m[max_j][max_i] = 0;
}

void jacobi_rotation2(double m[SIZE][SIZE]) {
  double max_nondiag = 0.0;
  int max_i = 0;
  int max_j = 0;
  for (int i = 0; i < SIZE; ++i)
    for (int j = i + 1; j < SIZE; ++j)
      if (fabs(m[i][j]) > max_nondiag) {
        max_i = i;
        max_j = j;
        max_nondiag = fabs(m[i][j]);
      }
  double const cur_ii = m[max_i][max_i];
  double const cur_jj = m[max_j][max_j];
  double const theta = 0.5 * atan(2 * m[max_i][max_j] / (cur_ii - cur_jj));
  double const c = cos(theta);
  double const s = sin(theta);
  double rotation_matrix[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      rotation_matrix[i][j] = (i == j) ? 1 : 0;
  rotation_matrix[max_i][max_i] =  c;
  rotation_matrix[max_j][max_j] =  c;
  rotation_matrix[max_i][max_j] = -s;
  rotation_matrix[max_j][max_i] =  s;
  printf("\nRotation matrix:\n");
  print_matrix(rotation_matrix);
  double new_matrix[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j) {
      new_matrix[i][j] = 0;
      for (int k = 0; k < SIZE; ++k)
        new_matrix[i][j] += m[i][k] * rotation_matrix[k][j];
    }
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j) {
      m[i][j] = 0;
      for (int k = 0; k < SIZE; ++k)
        m[i][j] += rotation_matrix[k][i] * new_matrix[k][j];
    }
}

double jacobi_norm(double m[SIZE][SIZE]) {
  double acc = 0;
  for (int i = 0; i < SIZE; ++i)
    for (int j = i + 1; j < SIZE; ++j)
      acc += fabs(m[i][j]);
  return acc;
}

void jacobi() {
  double jacobi_matrix[SIZE][SIZE] = { { 0 } };
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      jacobi_matrix[i][j] = matrix[i][j];
  for (int i = 1; jacobi_norm(jacobi_matrix) > EPS;  ++i) {
    jacobi_rotation(jacobi_matrix);
    printf("\n%d:\n", i);
    print_matrix(jacobi_matrix);
  }
}

void jacobi2() {
  double jacobi_matrix[SIZE][SIZE] = { { 0 } };
  for (int i = 0; i < SIZE; ++i)
    for (int j = 0; j < SIZE; ++j)
      jacobi_matrix[i][j] = matrix[i][j];
  for (int i = 1;;  ++i) {
    double diags[SIZE];
    for (int i = 0; i < SIZE; ++i)
      diags[i] = jacobi_matrix[i][i];
    jacobi_rotation2(jacobi_matrix);
    printf("\n%d:\n", i);
    print_matrix(jacobi_matrix);
    int done = 1;
    for (int i = 0; i < SIZE; ++i)
      if (fabs(diags[i] - jacobi_matrix[i][i]) > EPS) {
        done = 0;
        break;
      }
    if (done)
      break;
  }
}

int main() {
  printf("A = \n");
  print_matrix(matrix);
  printf("\ny = \n");
  print_vector(vector);
  printf("\n");

  int iter_counter = 0;
  double const max_lambda = iteration(&iter_counter);
  printf("Iteration: max by absolute value eigenvalue is %.*f (calculated in %d iterations)\n",
      DOUBLE_PREC, max_lambda, iter_counter);
  printf("\n");

  printf("\nJacobi:\n");
  jacobi();

  return 0;
}
