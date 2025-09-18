#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to evaluate polynomial f(x) = c0 + c1*x + c2*x^2 + ... + cN*x^N
double poly_eval(double *coeff, int N, double x) {
    double result = 0.0;
    for (int i = 0; i <= N; i++) {
        result += coeff[i] * pow(x, i);
    }
    return result;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <N> <c0 c1 ... cN> <M> <stddev>\n", argv[0]);
        return 1;
    }

    // Read N
    int N = atoi(argv[1]);

    // Read coefficients
    double *coeff = (double *)malloc((N + 1) * sizeof(double));
    for (int i = 0; i <= N; i++) {
        coeff[i] = atof(argv[2 + i]);
    }

    // Read M and stddev
    int M = atoi(argv[2 + N + 1]);
    double stddev = atof(argv[2 + N + 2]);

    // Generate 2M points in [0,1] (you can change interval if needed)
    int total_points = 2 * M;
    for (int i = 0; i < total_points; i++) {
        double x = (double)i / (total_points - 1); // evenly spaced in [0,1]
        double y_true = poly_eval(coeff, N, x);
        printf("%f %f\n", x, y_true);
    }

    free(coeff);
    return 0;
}

