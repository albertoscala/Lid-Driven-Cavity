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



int main() {
    cout << "Hello World!" << endl;
    return 0;
}