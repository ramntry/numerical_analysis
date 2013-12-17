#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "../support/vector.h"
#include "../support/polynomial.h"
#include "../hw03/finite_diff_table.h"

double pi;
double eps = 1.e-6;

int main(int argc, char **argv)
{
    srand(time(NULL));
    pi = 4.0 * atan(1);

    return 0;
}
