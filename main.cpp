//
//  advanced-programming
//
//  Created by Daniel Tompkins on 19/03/2025.
//

#include <iostream> // include iostream to allow inputs and outputs
#include "wave_solver.hpp"  // include solver class

// set up main function
int main() {
    // simulation parameters saved as double precision floats:
    // wave speed in km/s
    const double c = 4.0;
    // grid spacing in km (so spacing is 100m)
    const double dx = 0.1;
    const double dy = 0.1;
    // timestep in s
    const double dt = 0.01;
    //simulation parameters saved as integers
    // number of grid points - just a number, no unit attached
    const int nx = 100;
    const int ny = 100;
    // number of timesteps to run simulation for
    const int nt = 500;

    // create an instance of the solver class
    WaveSolver2D solver(c, dx, dy, dt, nx, ny, nt);
    // call function to initialise solution
    solver.initialiseSolution();
    // call function to solve for given parameters
    solver.solve();
    // call function to save the full solution data - must pass a path to the name/location you wish to save it is
    solver.saveSolution("/Users/danieltompkins/Documents/programming-coursework/solution.txt");
    // confirmation that the solution has been saved
    std::cout << "Solution saved to solution.txt" << std::endl;
    
    return 0;
}

