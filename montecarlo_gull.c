#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "random.c"
#define NUM_TRIES 10000000  // Number of Monte Carlo iterations
#define K 18822.79734      // Carrying capacity
#define ALPHA 0.3489494104672237     // Population growth rate without dispersal
#define BETA 2.43826356697e-5 // Intrinsic growth rate
#define EPSILON 0.11       // Death rate

// Structure to hold parameters
typedef struct {
    double x0;
    double phi;
    double lambda;
    double mu;
    double sigma;
    double delta;
} Parameters;

// Function to generate a random number in a given range
double rand_range(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}
// Function to compute the dispersal function phi(x)
double dispersal_function(double x, double mu, double sigma, double delta) {
    double E_x = sigma * (x - delta) / (1000 + sigma * fabs(x - delta));
    double E_dir = ((mu / (1000 + sigma * delta)) * (1 - x / delta) + (x / delta)) * E_x;
    return (x <= delta) ? (1 - E_dir) / (1 - E_dir) : (1 - E_x) / (1 - E_dir);
}

void dxdt(double t, double x, double *dx, void *params) {
    double *p = (double *)params;
    double phi = p[0];
    double lambda = p[1];
    double mu = p[2];
    double sigma = p[3];
    double delta = p[4];

    double dispersal = lambda * dispersal_function(x, mu, sigma, delta);
    *dx = phi * x - BETA * x * x - dispersal;
}



// Function to evaluate F(x0, phi, lambda, mu, sigma, delta)
double evaluate_F(Parameters p, double observed[], int years) {
    double x = p.x0;
    double error = 0.0;

    for (int t = 0; t < years; t++) {
        double dispersal = p.lambda * dispersal_function(x, p.mu, p.sigma, p.delta);
        x = (p.phi * x - BETA * x * x - dispersal);
        if (x < 0) x = 0; // Population cannot be negative
        error += pow(x - observed[t], 2);
    }
    return sqrt(error);
}

 


int main() {
    randomize();

    double observed[] = {15329, 14177, 13031, 9762, 11271, 8688, 7571, 6983, 4778, 2067, 1586, 793};
    int years = sizeof(observed) / sizeof(observed[0]);
    
    Parameters best_params;
    double best_F = INFINITY;
    double x0_step = 5206.0 / (1 << 19) - 1;
    double phi_step = (ALPHA - 0.12) / ((1ULL << 32) - 1);
    double lambda_step = 2700.0 / ((1 << 19) - 1);
    double mu_step = 10.0 / ((1 << 24) - 1);
    double sigma_step = 50.0 / ((1 << 19) - 1);
    long double delta_step = 20000.0 / ((1 << 16) - 1);
    // Monte Carlo Simulation
    for (int i = 0; i < NUM_TRIES; i++) {
        Parameters p = {
            discrete_random(19,12726,x0_step),  // x0
            discrete_random(32,0.12, phi_step),  // phi
            discrete_random(19,300, lambda_step),  // lambda
            discrete_random(24,0,  mu_step),  // mu
            discrete_random(19,0,  sigma_step),  // sigma
            discrete_random(16,0, delta_step)  // delta
        };
        
        double F_value = evaluate_F(p, observed, years);
        if (F_value < best_F) {
            best_F = F_value;
            best_params = p;
        }
    }

    // Print best found parameters
    printf("Best parameters found:\n");
    printf("x0 = %lf\n", best_params.x0);
    printf("phi = %lf\n", best_params.phi);
    printf("lambda = %lf\n", best_params.lambda);
    printf("mu = %lf\n", best_params.mu);
    printf("sigma = %lf\n", best_params.sigma);
    printf("delta = %lf\n", best_params.delta);
    printf("Best F value = %lf\n", best_F);

    return 0;
} 