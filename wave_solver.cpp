//
//  wave_solver.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 10/03/2025.
//

#include "wave_solver.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>  // For atof and atoi if needed

// Constructor implementation: initializes parameters and allocates memory.
WaveSolver2D::WaveSolver2D(double c, double dx, double dy, double dt, int nx, int ny, int nt)
    : c(c), dx(dx), dy(dy), dt(dt), nx(nx), ny(ny), nt(nt) {
    // Allocate memory for the solution: u[time][x][y]
    u.resize(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
    u_prev.resize(nx, std::vector<double>(ny, 0.0));
}

// Initialize the solution grid with a Gaussian pulse.
void WaveSolver2D::initialiseSolution() {
    double x0 = nx * dx / 2.0, y0 = ny * dy / 2.0; // Center of the pulse
    double sigma = 1.0;                             // Width of the pulse
    double A = 1.0;                                 // Amplitude

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = i * dx;
            double y = j * dy;
            u[0][i][j] = A * exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (sigma * sigma));
            u_prev[i][j] = u[0][i][j]; // Initial velocity is zero.
        }
    }
}

void WaveSolver2D::solve() {
    double C = c * dt / dx; // CFL number

    // Loop over time steps (starting from time index 1)
    for (int n = 1; n < nt - 1; ++n) {
        // Dynamically allocate a temporary 2D array for the next time step.
        double** u_next = new double*[nx];
        for (int i = 0; i < nx; ++i) {
            u_next[i] = new double[ny];
            // Initialize the row to zero (optional, but safe)
            for (int j = 0; j < ny; ++j)
                u_next[i][j] = 0.0;
        }

        // Compute the interior grid points using the wave equation update.
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                u_next[i][j] = 2 * u[n][i][j] - u[n - 1][i][j] +
                               C * C * (u[n][i + 1][j] + u[n][i - 1][j] +
                                        u[n][i][j + 1] + u[n][i][j - 1] - 4 * u[n][i][j]);
            }
        }

        // Apply homogeneous Dirichlet boundary conditions.
        for (int i = 0; i < nx; ++i) {
            u_next[i][0] = 0.0;          // Bottom boundary
            u_next[i][ny - 1] = 0.0;       // Top boundary
        }
        for (int j = 0; j < ny; ++j) {
            u_next[0][j] = 0.0;          // Left boundary
            u_next[nx - 1][j] = 0.0;       // Right boundary
        }

        // Copy the temporary array into the solution vector for the next time step.
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                u[n + 1][i][j] = u_next[i][j];
            }
        }

        // Deallocate the temporary 2D array.
        for (int i = 0; i < nx; ++i) {
            delete[] u_next[i];
        }
        delete[] u_next;
    }
}


// Save the solution to a file.
void WaveSolver2D::saveSolution(const std::string& filename) const {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    for (int n = 0; n < nt; ++n) {
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                outfile << u[n][i][j] << " ";
            }
            outfile << std::endl;
        }
        outfile << std::endl;
    }
    outfile.close();
}
