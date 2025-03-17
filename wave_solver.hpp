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

// the WaveSolver2D class is designed to solve the 2D wave equation using a finite-difference scheme
// it stores the solution on a 3D grid (time, x, y) and provides functions to initialise the solution,
// perform the time-stepping, and save the computed data to a file
class WaveSolver2D {
public:
    // constructor to initialise simulation parameters and allocate memory for the solution grids
    // parameters:
    //   c   - wave speed (km/s)
    //   dx  - grid spacing in the x-direction (km)
    //   dy  - grid spacing in the y-direction (km)
    //   dt  - time step (s)
    //   nx  - number of grid points in the x-direction
    //   ny  - number of grid points in the y-direction
    //   nt  - number of time steps to simulate
    WaveSolver2D(double c, double dx, double dy, double dt, int nx, int ny, int nt);

    // initialise the solution grid with a Gaussian pulse centred in the domain
    // this function sets the initial condition for the simulation
    void initialiseSolution();

    // solve the 2D wave equation using a finite-difference scheme
    // the solution is updated over time and stored in the 3D grid 'u'
    void solve();

    // save the complete solution data to a text file
    // the file is written in blocks, with each block corresponding to one time step
    void saveSolution(const std::string& filename) const;

private:
    // simulation parameters:
    double c;          // wave speed (km/s)
    double dx, dy;     // spatial grid spacing (km)
    double dt;         // time step (s)
    int nx, ny;        // number of grid points in x and y directions
    int nt;            // number of time steps

    // 3D grid for the solution data: indexed as [time][x][y]
    std::vector<std::vector<std::vector<double>>> u;

    // 2D grid for the solution at the previous time step, used in the time-stepping scheme
    std::vector<std::vector<double>> u_prev;
};

#endif // wave_solver_hpp


