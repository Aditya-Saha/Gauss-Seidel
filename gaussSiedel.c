#include <stdio.h>
#include <math.h>

#define N 3 // Size of the square matrix

// Function prototypes
void printMatrix(double A[N][N], double B[N]);
void gaussSeidel(double A[N][N], double B[N], double X[N], double tolerance, int maxIterations);

int main() {
    // Example usage:
    double A[N][N] = {{4, 1, 2}, {3, 5, 1}, {1, 1, 3}}; // Coefficients matrix
    double B[N] = {4, 7, 3}; // Constants vector
    double X[N] = {0}; // Initial guess for the solution
    double tolerance = 1e-6; // Tolerance for convergence
    int maxIterations = 1000; // Maximum number of iterations

    printf("System of equations:\n");
    printMatrix(A, B);

    gaussSeidel(A, B, X, tolerance, maxIterations);

    return 0;
}

// Function to print the matrix
void printMatrix(double A[N][N], double B[N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2lf ", A[i][j]);
        }
        printf("| %.2lf\n", B[i]);
    }
    printf("\n");
}

// Function to perform Gauss-Seidel iteration
void gaussSeidel(double A[N][N], double B[N], double X[N], double tolerance, int maxIterations) {
    double newX[N]; // Store updated values of X

    int iterations = 0;
    double error = tolerance + 1;

    while (error > tolerance && iterations < maxIterations) {
        error = 0.0;

        for (int i = 0; i < N; i++) {
            double sum = 0.0;

            for (int j = 0; j < N; j++) {
                if (j != i) {
                    sum += A[i][j] * X[j];
                }
            }

            newX[i] = (B[i] - sum) / A[i][i];

            error += fabs(newX[i] - X[i]);
            X[i] = newX[i];
        }

        iterations++;
    }

    if (iterations >= maxIterations) {
        printf("Maximum iterations reached. Solution may not be accurate.\n");
    } else {
        printf("Gauss-Seidel method converged in %d iterations.\n", iterations);
    }

    printf("Solution:\n");
    for (int i = 0; i < N; i++) {
        printf("X[%d] = %.6lf\n", i, X[i]);
    }
}
