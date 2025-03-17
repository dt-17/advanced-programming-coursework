//
//  wave_solver.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 10/03/2025.
//

#include "wave_solver.hpp"      // include header for the solver class
#include <cmath>                // include cmath for mathematical functions (e.g. exp)
#include <fstream>              // include fstream for file I/O
#include <iostream>             // include iostream for console I/O
#include <cstdlib>              // include cstdlib for atof and atoi if needed

// constructor for the WaveSolver2D class
// initialises simulation parameters and allocates memory for the solution arrays
WaveSolver2D::WaveSolver2D(double c, double dx, double dy, double dt, int nx, int ny, int nt)
    : c(c), dx(dx), dy(dy), dt(dt), nx(nx), ny(ny), nt(nt) {
    // allocate memory for the 3D solution vector: u[time][x][y]
    u.resize(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
    // allocate memory for the previous time step array (used in the time-stepping scheme)
    u_prev.resize(nx, std::vector<double>(ny, 0.0));
}

// initialise the solution grid with a Gaussian pulse centred in the domain
// this sets up the initial condition for the simulation
void WaveSolver2D::initialiseSolution() {
    double x0 = nx * dx / 2.0, y0 = ny * dy / 2.0; // calculate centre of the domain (pulse centre)
    double sigma = 1.0;                             // width of the Gaussian pulse
    double A = 1.0;                                 // amplitude of the pulse

    // loop over each grid point in the spatial domain
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double x = i * dx;
            double y = j * dy;
            // set the initial condition for time index 0 using the Gaussian function
            u[0][i][j] = A * exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0)) / (sigma * sigma));
            // copy the initial condition into the previous time step array (initial velocity is zero)
            u_prev[i][j] = u[0][i][j];
        }
    }
}

// solve the 2D wave equation using a finite difference scheme
// the solution is updated in time and stored in the 3D vector 'u'
void WaveSolver2D::solve() {
    double C = c * dt / dx; // calculate the CFL number - important part of timestepping scheme

    // loop over time steps (starting at time index 1 up to nt-2)
    for (int n = 1; n < nt - 1; ++n) {
        // dynamically allocate a temporary 2D array for the next time step
        double** u_next = new double*[nx];
        for (int i = 0; i < nx; ++i) {
            u_next[i] = new double[ny];
            // initialise each row to zero (ensures safety in computation)
            for (int j = 0; j < ny; ++j)
                u_next[i][j] = 0.0;
        }

        // compute the solution for the interior grid points using the update formula
        // the wave equation is discretised in both space and time
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                u_next[i][j] = 2 * u[n][i][j] - u[n - 1][i][j] +
                               C * C * (u[n][i + 1][j] + u[n][i - 1][j] +
                                        u[n][i][j + 1] + u[n][i][j - 1] - 4 * u[n][i][j]);
            }
        }

        // apply homogeneous Dirichlet boundary conditions
        // these conditions set the solution to zero along the boundaries
        for (int i = 0; i < nx; ++i) {
            u_next[i][0] = 0.0;          // bottom boundary
            u_next[i][ny - 1] = 0.0;       // top boundary
        }
        for (int j = 0; j < ny; ++j) {
            u_next[0][j] = 0.0;          // left boundary
            u_next[nx - 1][j] = 0.0;       // right boundary
        }

        // copy the computed values from the temporary array into the solution vector for the next time step
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                u[n + 1][i][j] = u_next[i][j];
            }
        }

        // deallocate the temporary 2D array to free memory
        for (int i = 0; i < nx; ++i) {
            delete[] u_next[i];
        }
        delete[] u_next;
    }
}

// save the complete solution data to a text file
// the data is written in a format with each time step's grid on a separate block
void WaveSolver2D::saveSolution(const std::string& filename) const {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    // loop over all time steps
    for (int n = 0; n < nt; ++n) {
        // loop over each row in the spatial domain
        for (int i = 0; i < nx; ++i) {
            // loop over each column in the row
            for (int j = 0; j < ny; ++j) {
                outfile << u[n][i][j] << " ";  // write each value followed by a space
            }
            outfile << std::endl;              // newline after each row
        }
        outfile << std::endl;                  // extra newline to separate time steps
    }
    outfile.close();                           // close the file
}
