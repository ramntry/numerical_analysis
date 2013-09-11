#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define DISPOSED_POLYNOMIAL_DEG (-1)
#define READ_POLYNOMIAL_ERROR (-1)
#define READ_POLYNOMIAL_OK (1)
#define PRINT_POLYNOMIAL_EPS (1.e-8)
#define PRINT_POLYNOMIAL_PRECISION (2)

#define NUMERIC_EPS (1.e-8)
#define PRINT_NUMERIC_PRECISION (9)
#define NUM_MONTE_CARLO_POINTS (200)
#define MAX_STARTING_POINT_DEVIATION (0.1)

#define max(a, b) ((a) < (b) ? (b) : (a))

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

typedef struct
{
    int deg;
    double *coeffs;
}
Polynomial;

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

double calcValue(Polynomial const *polynomial, double point)
{
    assert(polynomial->deg >= 0);
    double result = polynomial->coeffs[polynomial->deg];
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
    while (fabs(currentValue = func(funcArg, root)) >= NUMERIC_EPS) {
        root -= currentValue / deriv(derivArg, root);
    }
    return root;
}

void improveAccuracyOfRoots(NewtonFunction func, void *funcArg
                          , NewtonFunction deriv, void *derivArg, Vector *startingPoints)
{
    for (size_t i = 0; i < startingPoints->size; ++i) {
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

    Vector startingPoints = findStartingPoints(&p, negativeRootsLowerBound, positiveRootsUpperBound);
    printf("Starting points: ");
    printVector(&startingPoints, stdout);
    putchar('\n');

    Polynomial dp = getDerivative(&p);
    printf("Derivative: ");
    printPolynomial(&dp, stdout);
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
    printf("Root: %.*lf\n", PRINT_NUMERIC_PRECISION, root);
    printf("Value: %.*lf\n", PRINT_NUMERIC_PRECISION, complexFunction(NULL, root));

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
