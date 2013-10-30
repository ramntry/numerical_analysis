#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "polynomial.h"

void initPolynomial(Polynomial *polynomial, int deg)
{
    assert(deg >= DISPOSED_POLYNOMIAL_DEG);
    polynomial->deg = deg;
    if (deg == DISPOSED_POLYNOMIAL_DEG) {
        polynomial->coeffs = NULL;
        return;
    }
    polynomial->coeffs = (double *)malloc(sizeof(double) * (deg + 1));
}

void disposePolynomial(Polynomial *polynomial)
{
    free(polynomial->coeffs);
    polynomial->coeffs = NULL;
    polynomial->deg = DISPOSED_POLYNOMIAL_DEG;
}

int readPolynomial(Polynomial *polynomial, FILE *fd)
{
    if (fscanf(fd, "%d", &polynomial->deg) != 1 || polynomial->deg < 0) {
        polynomial->deg = DISPOSED_POLYNOMIAL_DEG;
        return READ_POLYNOMIAL_ERROR;
    }
    initPolynomial(polynomial, polynomial->deg);
    for (int i = polynomial->deg; i >= 0; --i) {
        if (fscanf(fd, "%lf", &polynomial->coeffs[i]) != 1) {
            disposePolynomial(polynomial);
            return READ_POLYNOMIAL_ERROR;
        }
    }
    return READ_POLYNOMIAL_OK;
}

void printPolynomial(Polynomial const *polynomial, FILE *fd)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    for (int i = polynomial->deg; i >= 0; --i) {
        if (i == polynomial->deg ) {
            if (i == 0 || fabs(polynomial->coeffs[i] - 1.) >= PRINT_POLYNOMIAL_EPS) {
                fprintf(fd, "%.*lf", PRINT_POLYNOMIAL_PRECISION, polynomial->coeffs[i]);
            }
        } else if (fabs(polynomial->coeffs[i]) >= PRINT_POLYNOMIAL_EPS) {
            if (polynomial->coeffs[i] >= 0) {
                fprintf(fd, " + ");
            } else {
                fprintf(fd, " - ");
            }
            if (fabs(polynomial->coeffs[i] - 1.) >= PRINT_POLYNOMIAL_EPS || i == 0) {
                fprintf(fd, "%.*lf", PRINT_POLYNOMIAL_PRECISION, fabs(polynomial->coeffs[i]));
            }
        } else {
            continue;
        }
        if (i >= 1) {
            fprintf(fd, "x");
        }
        if (i > 1) {
            fprintf(fd, "^%d", i);
        }
    }
}

Polynomial getPolynomialFromUser()
{
    printf("Enter polynomial's degree and coefficients separated by spaces: ");
    Polynomial polynomial;
    if (readPolynomial(&polynomial, stdin) != READ_POLYNOMIAL_OK) {
        fprintf(stderr, "Polynomial reading error!\n");
        exit(EXIT_FAILURE);
    }
    return polynomial;
}

void multiplyByScalar(Polynomial *polynomial, double scalar)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    for (int i = 0; i <= polynomial->deg; ++i) {
        polynomial->coeffs[i] *= scalar;
    }
}

void makePositiveLeadingCoeff(Polynomial *polynomial)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    if (polynomial->coeffs[polynomial->deg] < 0.0) {
        multiplyByScalar(polynomial, -1.0);
    }
}

