#include<iostream>

using namespace std;

// Physical parameters
const double L = 1.0; // Dimension of the cavity
const double rho = 1.0; // Density
const double mu = 0.1; // Viscosity
const double U = 1.0; // Lid velocity

// Grid parameters
const int Nx = 101; // Number of grid points in x-direction
const int Ny = 101; // Number of grid points in y-direction
const double dx = L/(Nx-1); // Grid spacing in x-direction
const double dy = L/(Ny-1); // Grid spacing in y-direction

// Time parameters
const double dt = 0.001; // Time step
const int num_steps = 1000; // Number of time steps

// Variable matrices
double u[Nx][Ny]; // x-velocity
double v[Nx][Ny]; // y-velocity
double p[Nx][Ny]; // Pressure

void initialize();
void apply_boundary_conditions();
void solve_navier_stokes();

int main() {
    initialize();

    // Main simulation loop
    for (int time_step = 0; time_step < 2; ++time_step) {
        apply_boundary_conditions();
        solve_navier_stokes();
    
        // Print out u, v, p
        cout << "Time step: " << time_step << endl;
        cout << "u:" << endl;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                cout << u[i][j] << " ";
            }
            cout << endl;
        }

        cout << "v:" << endl;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                cout << v[i][j] << " ";
            }
            cout << endl;
        }

        cout << "p:" << endl;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                cout << p[i][j] << " ";
            }
            cout << endl;
        }

    }

    return 0;
}

void initialize() {
    // Initialize u, v, p
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            u[i][j] = 1.0;
            v[i][j] = 1.0;
            p[i][j] = 1.0;
        }
    }
}

void apply_boundary_conditions() {
    // Apply boundary conditions to u, v, p

    // Left boundary
    for (int j = 0; j < Ny; ++j) {
        u[0][j] = 0.0;
        v[0][j] = 0.0;
        p[0][j] = p[1][j];
    }

    // Right boundary
    for (int j = 0; j < Ny; ++j) {
        u[Nx-1][j] = 0.0;
        v[Nx-1][j] = 0.0;
        p[Nx-1][j] = p[Nx-2][j];
    }

    // Bottom boundary
    for (int i = 0; i < Nx; ++i) {
        u[i][0] = 0.0;
        v[i][0] = 0.0;
        p[i][0] = p[i][1];
    }

    // Top boundary
    for (int i = 0; i < Nx; ++i) {
        u[i][Ny-1] = U;
        v[i][Ny-1] = 0.0;
        p[i][Ny-1] = p[i][Ny-2];
    }

}

void solve_navier_stokes() {
    // Solve Navier-Stokes equations for u, v, p

    // Temporary arrays to store intermediate velocities
    double u_temp[Nx][Ny];
    double v_temp[Nx][Ny];

    // Update u and v using the discretized Navier-Stokes equations
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u_temp[i][j] = u[i][j] + dt * (
                -u[i][j] * (u[i][j] - u[i-1][j]) / dx
                - v[i][j] * (u[i][j] - u[i][j-1]) / dy
                - (1.0 / rho) * (p[i+1][j] - p[i][j]) / dx
                + mu * (
                    (u[i+1][j] - 2.0 * u[i][j] + u[i-1][j]) / (dx * dx)
                    + (u[i][j+1] - 2.0 * u[i][j] + u[i][j-1]) / (dy * dy)
                )
            );

            v_temp[i][j] = v[i][j] + dt * (
                -u[i][j] * (v[i][j] - v[i-1][j]) / dx
                - v[i][j] * (v[i][j] - v[i][j-1]) / dy
                - (1.0 / rho) * (p[i][j+1] - p[i][j]) / dy
                + mu * (
                    (v[i+1][j] - 2.0 * v[i][j] + v[i-1][j]) / (dx * dx)
                    + (v[i][j+1] - 2.0 * v[i][j] + v[i][j-1]) / (dy * dy)
                )
            );
        }
    }

    // Solve the pressure Poisson equation
    // (Here, we use a simple iterative method, e.g., the Jacobi method)
    for (int iter = 0; iter < maxIterations; ++iter) {
        double p_temp[Nx][Ny];

        // Update pressure field
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                p_temp[i][j] = 0.25 * (
                    p[i+1][j] + p[i-1][j] + p[i][j+1] + p[i][j-1]
                    - dx * dy * rho / (2.0 * dt) * (
                        (u_temp[i+1][j] - u_temp[i-1][j]) / (2.0 * dx)
                        + (v_temp[i][j+1] - v_temp[i][j-1]) / (2.0 * dy)
                    )
                );
            }
        }

        // Update pressure field
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                p[i][j] = p_temp[i][j];
            }
        }
    }

    // Update velocities using the pressure gradient
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            u[i][j] = u_temp[i][j] - (dt / dx) * (p[i+1][j] - p[i][j]);
            v[i][j] = v_temp[i][j] - (dt / dy) * (p[i][j+1] - p[i][j]);
        }
    }

}