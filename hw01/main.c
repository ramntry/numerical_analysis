#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../support/polynomial.h"

#define NUMERIC_EPS (1.e-9)
#define PRINT_NUMERIC_PRECISION (12)

#define NUM_MONTE_CARLO_POINTS (200)

#ifndef max
#define max(a, b) ((a) < (b) ? (b) : (a))
#endif

typedef struct
{
    size_t capacity;
    size_t size;
    double *values;
}
Vector;

void initVector(Vector *vector, size_t capacity)
{
    vector->capacity = capacity;
    vector->size = 0;
    if (capacity > 0) {
        vector->values = (double *)malloc(sizeof(double) * capacity);
    } else {
        vector->values = NULL;
    }
}

void disposeVector(Vector *vector)
{
    free(vector->values);
    vector->values = NULL;
    vector->size = 0;
    vector->capacity = 0;
}

void append(Vector *vector, double value)
{
    assert(vector->size < vector->capacity);
    vector->values[vector->size++] = value;
}

void printVector(Vector const *vector, FILE *fd)
{
    fprintf(fd, "[");
    if (vector->size != 0) {
        fprintf(fd, "%.*lf", PRINT_NUMERIC_PRECISION, vector->values[0]);
    }
    for (size_t i = 1; i < vector->size; ++i) {
        fprintf(fd, ", %.*lf", PRINT_NUMERIC_PRECISION, vector->values[i]);
    }
    fprintf(fd, "]");
}

double rootsUpperBound(Polynomial const *polynomial)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    assert(polynomial->coeffs != NULL && polynomial->coeffs[polynomial->deg] > 0.0);
    double maxAbsOfNegativeCoeff = 0.0;
    for (int i = 0; i <= polynomial->deg; ++i) {
        if (polynomial->coeffs[i] < 0.0) {
            maxAbsOfNegativeCoeff = max(maxAbsOfNegativeCoeff, -polynomial->coeffs[i]);
        }
    }
    if (maxAbsOfNegativeCoeff == 0.0) {
        return NUMERIC_EPS;
    }
    int numOfFirstNegativeCoeff = 1;
    for (;; ++numOfFirstNegativeCoeff) {
        if (polynomial->coeffs[polynomial->deg - numOfFirstNegativeCoeff] < 0.0) {
            break;
        }
    }
    return 1.0 + pow(maxAbsOfNegativeCoeff / polynomial->coeffs[polynomial->deg], 1.0 / numOfFirstNegativeCoeff);
}

void createReverse(Polynomial const *src, Polynomial *dst)
{
    assert(src->deg != DISPOSED_POLYNOMIAL_DEG && src->deg == dst->deg);
    assert(src->coeffs != NULL && dst->coeffs != NULL);
    for (int i = 0; i <= src->deg; ++i) {
        dst->coeffs[i] = src->coeffs[src->deg - i];
    }
}

void createReflected(Polynomial const *src, Polynomial *dst)
{
    assert(src->deg != DISPOSED_POLYNOMIAL_DEG && src->deg == dst->deg);
    assert(src->coeffs != NULL && dst->coeffs != NULL);
    for (int i = 0; i <= src->deg; ++i) {
        dst->coeffs[i] = ((i % 2 == 0) ? 1.0 : -1.0) * src->coeffs[i];
    }
}

double rootsLowerBound(Polynomial const *polynomial)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    assert(polynomial->coeffs != NULL);
    assert(polynomial->coeffs[polynomial->deg] > 0.0);
    Polynomial reflected;
    initPolynomial(&reflected, polynomial->deg);
    createReflected(polynomial, &reflected);
    makePositiveLeadingCoeff(&reflected);
    double lowerBound = -rootsUpperBound(&reflected);
    disposePolynomial(&reflected);
    return lowerBound;
}

Vector findStartingPoints(Polynomial const *polynomial, double lowerBound, double upperBound)
{
    double const step = (upperBound - lowerBound) / NUM_MONTE_CARLO_POINTS;
    Vector startingPoints;
    initVector(&startingPoints, polynomial->deg);
    double prevValue = calcValue(polynomial, lowerBound);
    for (double point = lowerBound + step; point <= upperBound; point += step) {
        double const currentValue = calcValue(polynomial, point);
        if (prevValue * currentValue <= 0.0) {
            append(&startingPoints, point);
        }
        prevValue = currentValue;
    }
    return startingPoints;
}

typedef double (*NewtonFunction)(void *additionalArg, double point);

double newtonRoot(NewtonFunction func, void *funcArg, NewtonFunction deriv, void *derivArg, double startingPoint)
{
    double root = startingPoint;
    double currentValue = 0.0;
    int numOfStep = 0;
    while (fabs(currentValue = func(funcArg, root)) >= NUMERIC_EPS) {
        printf("step: %-5dx: %-20.*lff(x): %-20.*lf\n", ++numOfStep, PRINT_NUMERIC_PRECISION, root
                                                      , PRINT_NUMERIC_PRECISION, currentValue);
        root -= currentValue / deriv(derivArg, root);
    }
    printf("step: %-5dx: %-20.*lff(x): %-20.*lf\n", ++numOfStep, PRINT_NUMERIC_PRECISION, root
                                                  , PRINT_NUMERIC_PRECISION, func(funcArg, root));
    return root;
}

void improveAccuracyOfRoots(NewtonFunction func, void *funcArg
                          , NewtonFunction deriv, void *derivArg, Vector *startingPoints)
{
    for (size_t i = 0; i < startingPoints->size; ++i) {
        if (i > 0) {
            putchar('\n');
        }
        startingPoints->values[i] = newtonRoot(func, funcArg, deriv, derivArg, startingPoints->values[i]);
    }
}

double unsafeCalcValue(void *polynomial, double point)
{
    return calcValue((Polynomial *)polynomial, point);
}

int polynomialMain()
{
    Polynomial p = getPolynomialFromUser();
    printf("Your polynomial: ");
    printPolynomial(&p, stdout);
    putchar('\n');

    makePositiveLeadingCoeff(&p);
    double const negativeRootsLowerBound = rootsLowerBound(&p);
    double const positiveRootsUpperBound = rootsUpperBound(&p);
    printf("Lower bound of negative roots: %lf\n", negativeRootsLowerBound);
    printf("Upper bound of positive roots: %lf\n", positiveRootsUpperBound);

    Polynomial dp = getDerivative(&p);
    printf("Derivative: ");
    printPolynomial(&dp, stdout);
    putchar('\n');

    Vector startingPoints = findStartingPoints(&p, negativeRootsLowerBound, positiveRootsUpperBound);
    printf("Starting points: ");
    printVector(&startingPoints, stdout);
    putchar('\n');

    improveAccuracyOfRoots(unsafeCalcValue, &p, unsafeCalcValue, &dp, &startingPoints);
    printf("Roots: ");
    printVector(&startingPoints, stdout);
    putchar('\n');

    Vector values = calcValues(&p, &startingPoints);
    printf("Values: ");
    printVector(&values, stdout);
    putchar('\n');

    disposeVector(&values);
    disposeVector(&startingPoints);
    disposePolynomial(&dp);
    disposePolynomial(&p);

    return 0;
}

double scanStartingPoint(int argc, char **argv)
{
    assert(argc > 1);
    double startingPoint = 0.0;
    if (sscanf(argv[1], "%lf", &startingPoint) != 1) {
        fprintf(stderr, "Starting point scaning error (bad command line argument)!\n");
        exit(EXIT_FAILURE);
    }
    return startingPoint;
}

double complexFunction(void *additionalArg, double x)
{
    return x*x + 4*sin(x) - 1;
}

double derivOfComplexFunction(void *additionalArg, double x)
{
    return 2*x + 4*cos(x);
}

int complexFunctionMain(int argc, char **argv)
{
    double const startingPoint = scanStartingPoint(argc, argv);
    printf("Your starting point: %lf\n", startingPoint);

    double const root = newtonRoot(complexFunction, NULL, derivOfComplexFunction, NULL, startingPoint);
    printf("Roots: %.*lf\n", PRINT_NUMERIC_PRECISION, root);
    printf("Values: %.*lf\n", PRINT_NUMERIC_PRECISION, complexFunction(NULL, root));

    return 0;
}

int isSmallDegree(Polynomial const *polynomial)
{
    assert(polynomial->deg >= 0);
    return polynomial->deg < 3;
}

Vector smallDegreeSolve(Polynomial *const polynomial)
{
    assert(0 <= polynomial->deg && polynomial->deg < 3);
    assert(polynomial->coeffs != NULL && polynomial->coeffs[polynomial->deg] > 0.0);
    Vector roots;
    if (polynomial->deg == 0) {
        initVector(&roots, 0);
        return roots;
    }
    if (polynomial->deg == 1) {
        initVector(&roots, 1);
        append(&roots, -polynomial->coeffs[0] / polynomial->coeffs[1]);
        return roots;
    }
    double const c = polynomial->coeffs[0];
    double const b = polynomial->coeffs[1];
    double const a = polynomial->coeffs[2];
    double const discriminant = b*b - 4*a*c;
    if (discriminant < -NUMERIC_EPS/2) {
        initVector(&roots, 0);
        return roots;
    }
    double const firstPart = -0.5 * b / a;
    if (discriminant >= NUMERIC_EPS/2) {
        double const secondPart = 0.5 * sqrt(discriminant) / a;
        initVector(&roots, 2);
        append(&roots, firstPart - secondPart);
        append(&roots, firstPart + secondPart);
        return roots;
    }
    initVector(&roots, 1);
    append(&roots, firstPart);
    return roots;
}

Vector recursiveSolve(Polynomial *const polynomial)
{
    assert(polynomial->deg != DISPOSED_POLYNOMIAL_DEG);
    assert(polynomial->coeffs != NULL);
    assert(polynomial->coeffs[polynomial->deg] > 0.0);

    if (isSmallDegree(polynomial)) {
        return smallDegreeSolve(polynomial);
    }
    Polynomial derivative = getDerivative(polynomial);
    Vector rootsOfDeriv = recursiveSolve(&derivative);
    Vector roots;
    initVector(&roots, rootsOfDeriv.size + 1);

    double const upperBound = rootsUpperBound(polynomial);
    double prevPoint = rootsLowerBound(polynomial);
    double prevValue = calcValue(polynomial, prevPoint);
    for (int i = 0; i < rootsOfDeriv.size; ++i) {
        double const currentValue = calcValue(polynomial, rootsOfDeriv.values[i]);
        if (prevValue * currentValue <= 0.0) {
            append(&roots, 0.5 * (rootsOfDeriv.values[i] + prevPoint));
        }
        prevPoint = rootsOfDeriv.values[i];
        prevValue = currentValue;
    }
    if (prevValue * calcValue(polynomial, upperBound) <= 0.0) {
        append(&roots, 0.5 * (rootsOfDeriv.values[rootsOfDeriv.size - 1] + upperBound));
    }

    improveAccuracyOfRoots(unsafeCalcValue, polynomial, unsafeCalcValue, &derivative, &roots);

    disposeVector(&rootsOfDeriv);
    disposePolynomial(&derivative);
    return roots;
}

int recursiveMain()
{
    Polynomial p = getPolynomialFromUser();
    Vector roots = recursiveSolve(&p);

    printf("Your polynomial: ");
    printPolynomial(&p, stdout);
    printf("\nRoots: ");
    printVector(&roots, stdout);
    putchar('\n');

    disposeVector(&roots);
    disposePolynomial(&p);

    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        return polynomialMain();
    }
    if (strcmp(argv[1], "-r") == 0) {
        return recursiveMain();
    }
    return complexFunctionMain(argc, argv);
}
