//
//  wave_solver.hpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 10/03/2025.
//

#ifndef wave_solver_hpp
#define wave_solver_hpp

#include <vector>
#include <string>

class WaveSolver2D {
public:
    // Constructor to initialize parameters
    WaveSolver2D(double c, double dx, double dy, double dt, int nx, int ny, int nt);

    // Initialize the solution grid with a Gaussian pulse
    void initializeSolution();

    // Solve the 2D wave equation using a finite-difference scheme
    void solve();

    // Save the complete solution to a file
    void saveSolution(const std::string& filename) const;

private:
    // Member variables
    double c;          // Wave speed
    double dx, dy;     // Spatial steps
    double dt;         // Time step
    int nx, ny;        // Number of grid points
    int nt;            // Number of time steps

    // 3D grid for the solution: time, x, y
    std::vector<std::vector<std::vector<double>>> u;

    // 2D grid for the solution at the previous time step
    std::vector<std::vector<double>> u_prev;
};

#endif // WAVESOLVER2D_H

