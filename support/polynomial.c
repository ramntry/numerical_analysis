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
    polynomial->coeffs = (long double *)malloc(sizeof(long double) * (deg + 1));
    polynomial->effective_deg = UNKNOWN_EFFECTIVE_DEG;
}

void disposePolynomial(Polynomial *polynomial)
{
    free(polynomial->coeffs);
    polynomial->coeffs = NULL;
    polynomial->deg = DISPOSED_POLYNOMIAL_DEG;
    polynomial->effective_deg = UNKNOWN_EFFECTIVE_DEG;
}

int readPolynomial(Polynomial *polynomial, FILE *fd)
{
    if (fscanf(fd, "%d", &polynomial->deg) != 1 || polynomial->deg < 0) {
        polynomial->deg = DISPOSED_POLYNOMIAL_DEG;
        return READ_POLYNOMIAL_ERROR;
    }
    initPolynomial(polynomial, polynomial->deg);
    for (int i = polynomial->deg; i >= 0; --i) {
        if (fscanf(fd, "%Lf", &polynomial->coeffs[i]) != 1) {
            disposePolynomial(polynomial);
            return READ_POLYNOMIAL_ERROR;
        }
    }
    polynomial->effective_deg = polynomial->deg;
    return READ_POLYNOMIAL_OK;
}

void printPolynomial(Polynomial const *polynomial, FILE *fd)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    for (int i = polynomial->deg; i >= 0; --i) {
        if (i == polynomial->deg ) {
            if (i == 0 || fabsl(polynomial->coeffs[i] - 1.) >= PRINT_POLYNOMIAL_EPS) {
                fprintf(fd, "%.*Lf", PRINT_POLYNOMIAL_PRECISION, polynomial->coeffs[i]);
            }
        } else if (fabsl(polynomial->coeffs[i]) >= PRINT_POLYNOMIAL_EPS) {
            if (polynomial->coeffs[i] >= 0) {
                fprintf(fd, " + ");
            } else {
                fprintf(fd, " - ");
            }
            if (fabsl(polynomial->coeffs[i] - 1.) >= PRINT_POLYNOMIAL_EPS || i == 0) {
                fprintf(fd, "%.*Lf", PRINT_POLYNOMIAL_PRECISION, fabsl(polynomial->coeffs[i]));
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

void multiplyByScalar(Polynomial *polynomial, long double scalar)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    for (int i = 0; i <= polynomial->deg; ++i) {
        polynomial->coeffs[i] *= scalar;
    }
}

void setToScalar(Polynomial *polynomial, long double scalar)
{
    assert(polynomial->deg >= 0);
    polynomial->coeffs[0] = scalar;
    for (int i = 1; i <= polynomial->deg; ++i) {
        polynomial->coeffs[i] = 0.0;
    }
    polynomial->effective_deg = 0;
}

void addRoot(Polynomial *polynomial, long double root)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG
        && polynomial->effective_deg != UNKNOWN_EFFECTIVE_DEG
        && polynomial->effective_deg < polynomial->deg);
    int new_deg = ++polynomial->effective_deg;
    polynomial->coeffs[new_deg] = polynomial->coeffs[new_deg - 1];
    for (int i = new_deg - 1; i > 0; --i) {
        polynomial->coeffs[i] *= -root;
        polynomial->coeffs[i] += polynomial->coeffs[i - 1];
    }
    polynomial->coeffs[0] *= -root;
}

void addPolynomial(Polynomial *dst, Polynomial const *src)
{
    assert(dst->deg != DISPOSED_POLYNOMIAL_DEG
        && dst->effective_deg != UNKNOWN_EFFECTIVE_DEG
        && dst->effective_deg <= dst->deg
        && src->deg != DISPOSED_POLYNOMIAL_DEG
        && src->effective_deg != UNKNOWN_EFFECTIVE_DEG
        && src->effective_deg <= src->deg
        && dst->deg <= src->effective_deg);
    for (int i = dst->effective_deg + 1; i <= src->effective_deg; ++i) {
        dst->coeffs[i] = 0.0;
    }
    dst->effective_deg = src->effective_deg > dst->effective_deg
                       ? src->effective_deg : dst->effective_deg;
    for (int i = 0; i <= src->effective_deg; ++i) {
        dst->coeffs[i] += src->coeffs[i];
    }
}

void makePositiveLeadingCoeff(Polynomial *polynomial)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    if (polynomial->coeffs[polynomial->deg] < 0.0) {
        multiplyByScalar(polynomial, -1.0);
    }
}

Polynomial getDerivative(Polynomial const *polynomial)
{
    assert(polynomial->deg >= 0);
    Polynomial derivative;
    if (polynomial->deg == 0) {
        initPolynomial(&derivative, 0);
        derivative.coeffs[0] = 0.0;
        return derivative;
    }
    initPolynomial(&derivative, polynomial->deg - 1);
    for (int i = 1; i <= polynomial->deg; ++i) {
        derivative.coeffs[i - 1] = i * polynomial->coeffs[i];
    }
    return derivative;
}

long double calcValue(Polynomial const *polynomial, long double point)
{
    assert(polynomial->deg >= 0);
    long double result = polynomial->coeffs[polynomial->deg];
    for (int i = polynomial->deg - 1; i >= 0; --i) {
        result *= point;
        result += polynomial->coeffs[i];
    }
    return result;
}

Vector calcValues(Polynomial const *polynomial, Vector const *points)
{
    Vector values;
    initVector(&values, points->size);
    values.size = values.capacity;
    for (size_t i = 0; i < values.size; ++i) {
        values.values[i] = calcValue(polynomial, points->values[i]);
    }
    return values;
}

