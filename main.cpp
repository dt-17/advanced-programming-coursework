//
//  main.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 20/01/2025.
//

#include <iostream>
#include "solution_visualiser.hpp"
#include "wave_solver.hpp"  // Your solver class

int main() {
    // Simulation parameters.
    const double c = 3.0;
    const double dx = 0.1;
    const double dy = 0.1;
    const double dt = 0.01;
    const int nx = 100;
    const int ny = 100;
    const int nt = 500;

    // Create and run the solver.
    WaveSolver2D solver(c, dx, dy, dt, nx, ny, nt);
    solver.initializeSolution();
    solver.solve();
    solver.saveSolution("/Users/danieltompkins/Documents/programming-coursework/solution.txt");
    std::cout << "Solution saved to solution.txt" << std::endl;

    // Create a Visualizer instance.
    Visualizer viz("/Users/danieltompkins/Documents/programming-coursework/solution.txt", nx, ny, nt);
    if (!viz.loadSolution()) {
        std::cerr << "Error loading solution file." << std::endl;
        return -1;
    }
    
    // Save the output of the final time step solution
    viz.saveFinalTimeStep("/Users/danieltompkins/Documents/programming-coursework/final_timestep.png");
    
    // Animate the solution.
    // Assuming you want the domain to be from 0 to 10 in both x and y,
    // and a pause of 0.05 seconds between frames.
    // viz.showAnimation(0.0, 10.0, 0.0, 10.0, 0.01);

    return 0;
}
