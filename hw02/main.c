#include <assert.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

typedef double (*Function2D)(void const *additionalArg, double x, double y);

double const eps = 1.e-9;
double const startPointDistance = 0.1;

inline double max(double x, double y)
{
    return x < y ? y : x;
}

typedef struct
{
    double a;
    double k;
}
Parameters;

typedef struct
{
    double x;
    double y;
}
Point2D;

void printPoint(Point2D point)
{
    printf("(%.10lf, %.10lf)", point.x, point.y);
}

typedef struct
{
    Point2D leftBottom;
    Point2D rightTop;
}
Rectangle;

inline double square(double x)
{
    return x*x;
}

double firstFunction(void const *parameters, double x, double y)
{
    double const k = ((Parameters const *)parameters)->k;
    return tanh(square(x) - y) - k*(x + y);
}

double first_x(void const *parameters, double x, double y)
{
    double const k = ((Parameters const *)parameters)->k;
    return 2.0*x/square(cosh(square(x) - y)) - k;
}

double first_y(void const *parameters, double x, double y)
{
    double const k = ((Parameters const *)parameters)->k;
    return -1.0/square(cosh(square(x) - y)) - k;
}

double secondFunction(void const *parameters, double x, double y)
{
    double const a = ((Parameters const *)parameters)->a;
    return square(x - 0.2) - a*square(y) - 1.5;
}

double second_x(void const *parameters, double x, double y)
{
    return 2.0*(x - 0.2);
}

double second_y(void const *parameters, double x, double y)
{
    double const a = ((Parameters const *)parameters)->a;
    return -2.0*a*y;
}

Point2D findOneStartPoint(Function2D f, Function2D g, Parameters const *parameters
                        , Rectangle *where, unsigned granularity, double targetDistance)
{
    assert(granularity > 0);
    double const stepX = (where->rightTop.x - where->leftBottom.x) / (granularity + 1);
    double const stepY = (where->rightTop.y - where->leftBottom.y) / (granularity + 1);
    Point2D currentPoint = { where->leftBottom.x, where->leftBottom.y };
    Point2D bestPoint = { NAN, NAN };
    double bestDistance = DBL_MAX;
    for (unsigned i = 1; i <= granularity; ++i) {
        currentPoint.y += stepY;
        for (unsigned j = 1; j <= granularity; ++j) {
            currentPoint.x += stepX;
            double const distance = max(fabs(f(parameters, currentPoint.x, currentPoint.y))
                                      , fabs(g(parameters, currentPoint.x, currentPoint.y)));
            if (distance < bestDistance) {
                bestPoint = currentPoint;
                bestDistance = distance;
            }
        }
        currentPoint.x = where->leftBottom.x;
    }
    if (bestDistance > targetDistance) {
        where->leftBottom.x = bestPoint.x - stepX;
        where->leftBottom.y = bestPoint.y - stepY;
        where->rightTop.x = bestPoint.x + stepX;
        where->rightTop.y = bestPoint.y + stepY;
        return findOneStartPoint(f, g, parameters, where, granularity, targetDistance);
    }
    return bestPoint;
}

Point2D newton2D(Function2D f, Function2D f_x, Function2D f_y
               , Function2D g, Function2D g_x, Function2D g_y
               , Parameters const *parameters, Point2D const *startPoint, double eps)
{
    printf("a = %lf, k = %lf, start point = ", parameters->a, parameters->k);
    printPoint(*startPoint);

    double x = startPoint->x;
    double y = startPoint->y;
    double distance = DBL_MAX;
    int stepCounter = 0;
    while ((distance = max(fabs(f(parameters, x, y)), fabs(g(parameters, x, y)))) > eps) {
        double const fv = f(parameters, x, y);
        double const fx = f_x(parameters, x, y);
        double const fy = f_y(parameters, x, y);
        double const gv = g(parameters, x, y);
        double const gx = g_x(parameters, x, y);
        double const gy = g_y(parameters, x, y);
        double const jacobian = fx*gy - fy*gx;
        printf("\nstep:\t%d | x: %13.10lf | y: %13.10lf | f(x, y): %13.10lf | g(x, y): %13.10lf", ++stepCounter, x, y, fv, gv);

        x -= (fv*gy - fy*gv) / jacobian;
        y -= (fx*gv - fv*gx) / jacobian;
    }
    printf("\nstep:\t%d | x: %13.10lf | y: %13.10lf | f(x, y): %13.10lf | g(x, y): %13.10lf\n\n"
         , ++stepCounter, x, y, f(parameters, x, y), g(parameters, x, y));
    Point2D result = { x, y };
    return result;
}

void findAndPrintSolution(double a, double k)
{
    Parameters parameters = { a, k };
    Rectangle rectangle = { { 0.0, 0.0 }, { 2.0, 2.0 } };
    Point2D startPoint = findOneStartPoint(firstFunction, secondFunction, &parameters, &rectangle, 9, startPointDistance);

    newton2D(firstFunction, first_x, first_y, secondFunction, second_x, second_y, &parameters, &startPoint, eps);
}

int mainSolveAllAtSmallRect()
{
    int const numofA = 6;
    int const numofK = 5;
    double as[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    double ks[] = { 0.44, 0.46, 0.48, 0.50, 0.52 };
    for (int i = 0; i < numofK; ++i) {
        for (int j = 0; j < numofA; ++j) {
            printf("Equation %d: ", i * numofA + j + 1);
            findAndPrintSolution(as[j], ks[i]);
        }
    }
    return 0;
}

int main()
{
    return mainSolveAllAtSmallRect();
}
