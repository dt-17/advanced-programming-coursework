//
//  main.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 20/01/2025.
//

#include <iostream>
#include "wave_solver.hpp"  // solver class
#include "solution_visualiser.hpp" // visualisation class

int main() {
    // Simulation parameters:
    // wave speed in km/s
    const double c = 4.0;
    // grid spacing in km (so spacing is 100m)
    const double dx = 0.1;
    const double dy = 0.1;
    // timestep in s
    const double dt = 0.01;
    // number of grid points - just a number, no unit attached
    const int nx = 100;
    const int ny = 100;
    // number of timesteps to run simulation for
    const int nt = 500;

    // Create and run the solver.
    WaveSolver2D solver(c, dx, dy, dt, nx, ny, nt);
    // call function to initialise solution
    solver.initialiseSolution();
    // call function to solve for given parameters
    solver.solve();
    // call function to save the full solution data - must pass a path to the name/location you wish to save it is
    solver.saveSolution("/Users/danieltompkins/Documents/programming-coursework/solution1.txt");
    // confirmation that the solution has been saved
    std::cout << "Solution saved to solution.txt" << std::endl;

    // Create a Visualizer instance.
    Visualiser vis("/Users/danieltompkins/Documents/programming-coursework/solution1.txt", nx, ny, nt);
    if (!vis.loadSolution()) {
        std::cerr << "Error loading solution file." << std::endl;
        return -1;
    }
    
    // Save the output of the final time step solution
    vis.saveFinalTimeStep("/Users/danieltompkins/Documents/programming-coursework/final_timestep1.png");
    
    // Animate the solution.
    // Assuming you want the domain to be from 0 to 10 in both x and y,
    // and a pause of 0.05 seconds between frames.
    //viz.showAnimation(0.0, 10.0, 0.0, 10.0, 0.01);

    return 0;
}
