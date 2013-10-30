#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <stdio.h>

#define DISPOSED_POLYNOMIAL_DEG (-1)
#define UNKNOWN_EFFECTIVE_DEG (-2)
#define READ_POLYNOMIAL_ERROR (-1)
#define READ_POLYNOMIAL_OK (1)
#define PRINT_POLYNOMIAL_EPS (1.e-8)
#define PRINT_POLYNOMIAL_PRECISION (2)

typedef struct
{
    int deg;
    double *coeffs;
    int effective_deg;
}
Polynomial;

void initPolynomial(Polynomial *polynomial, int deg);
void disposePolynomial(Polynomial *polynomial);
int readPolynomial(Polynomial *polynomial, FILE *fd);
Polynomial getPolynomialFromUser();
void printPolynomial(Polynomial const *polynomial, FILE *fd);
void setToScalar(Polynomial *polynomial, double scalar);
void multiplyByScalar(Polynomial *polynomial, double scalar);
void makePositiveLeadingCoeff(Polynomial *polynomial);
void addRoot(Polynomial *polynomial, double root);

#endif
