#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Evaluate polynomial f(x) = a0 + a1*x + ... + aN*x^N
double evaluate_polynomial(double x, int N, double *coeffs) {
    double result = 0.0;
    for (int i = 0; i <= N; i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Generate Gaussian noise using Box-Muller transform
double gaussian_noise(double stddev) {
    double u1 = ((double) rand() + 1) / ((double) RAND_MAX + 2);
    double u2 = ((double) rand() + 1) / ((double) RAND_MAX + 2);
    return stddev * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

// Lagrange interpolation
double lagrange_interpolate(double *x, double *y, int n, double xp) {
    double yp = 0.0;
    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (j != i)
                term *= (xp - x[j]) / (x[i] - x[j]);
        }
        yp += term;
    }
    return yp;
}

// Newton interpolation
double newton_interpolate(double *x, double *y, int n, double xp) {
    double *div_diff = malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) div_diff[i] = y[i];

    for (int i = 1; i < n; i++) {
        for (int j = n - 1; j >= i; j--) {
            div_diff[j] = (div_diff[j] - div_diff[j - 1]) / (x[j] - x[j - i]);
        }
    }

    double yp = div_diff[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        yp = yp * (xp - x[i]) + div_diff[i];
    }

    free(div_diff);
    return yp;
}

// RMSE calculation
double compute_rmse(double *predicted, double *actual, int n) {
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = predicted[i] - actual[i];
        sum_sq += diff * diff;
    }
    return sqrt(sum_sq / n);
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        printf("Usage: %s N a0 a1 ... aN M stddev\n", argv[0]);
        return 1;
    }

    srand((unsigned int) time(NULL));

    int N = atoi(argv[1]);
    double *coeffs = malloc((N + 1) * sizeof(double));
    for (int i = 0; i <= N; i++) {
        coeffs[i] = atof(argv[2 + i]);
    }

    int M = atoi(argv[2 + N + 1]);
    double stddev = atof(argv[2 + N + 2]);
    int total_points = 2 * M;

    double *x = malloc(total_points * sizeof(double));
    double *y = malloc(total_points * sizeof(double));
    double *y_true = malloc(total_points * sizeof(double));

    // Generate noisy and noise-free data
    for (int i = 0; i < total_points; i++) {
        x[i] = (double) i / (total_points - 1);
        y_true[i] = evaluate_polynomial(x[i], N, coeffs);
        y[i] = y_true[i] + gaussian_noise(stddev);
    }

    // Allocate interpolation arrays
    double *x_interp = malloc(M * sizeof(double));
    double *y_interp = malloc(M * sizeof(double));
    double *x_skipped = malloc(M * sizeof(double));
    double *y_skipped_true = malloc(M * sizeof(double));
    double *y_skipped_lagrange = malloc(M * sizeof(double));
    double *y_skipped_newton = malloc(M * sizeof(double));

    // Split alternate and skipped points
    for (int i = 0; i < M; i++) {
        x_interp[i] = x[2 * i];
        y_interp[i] = y[2 * i];
        x_skipped[i] = x[2 * i + 1];
        y_skipped_true[i] = y_true[2 * i + 1];
    }

    // Interpolate skipped points
    for (int i = 0; i < M; i++) {
        y_skipped_lagrange[i] = lagrange_interpolate(x_interp, y_interp, M, x_skipped[i]);
        y_skipped_newton[i] = newton_interpolate(x_interp, y_interp, M, x_skipped[i]);
    }

    // Compute RMSE
    double rmse_lagrange = compute_rmse(y_skipped_lagrange, y_skipped_true, M);
    double rmse_newton   = compute_rmse(y_skipped_newton, y_skipped_true, M);

    printf("RMSE (Lagrange): %lf\n", rmse_lagrange);
    printf("RMSE (Newton):   %lf\n", rmse_newton);

    // Cleanup
    free(coeffs);
    free(x); free(y); free(y_true);
    free(x_interp); free(y_interp);
    free(x_skipped); free(y_skipped_true);
    free(y_skipped_lagrange); free(y_skipped_newton);

    return 0;
}

