#include <iostream>
#include "wave_solver.hpp"  // include solver class

int main() {
    // Declare variables for simulation parameters.
    double c, dx, dy, dt;
    int nx, ny, nt;

    // Prompt the user for simulation parameters.
    std::cout << "Enter wave speed (km/s): ";
    std::cin >> c;

    std::cout << "Enter grid spacing in x (km): ";
    std::cin >> dx;

    std::cout << "Enter grid spacing in y (km): ";
    std::cin >> dy;

    std::cout << "Enter time step (s): ";
    std::cin >> dt;

    std::cout << "Enter number of grid points in x: ";
    std::cin >> nx;

    std::cout << "Enter number of grid points in y: ";
    std::cin >> ny;

    std::cout << "Enter number of time steps: ";
    std::cin >> nt;

    // Create an instance of the WaveSolver2D with user input parameters.
    WaveSolver2D solver(c, dx, dy, dt, nx, ny, nt);

    // Initialize the solution with a Gaussian pulse.
    solver.initialiseSolution();

    // Solve the 2D wave equation.
    solver.solve();

    // Save the solution to file.
    solver.saveSolution("/Users/danieltompkins/Documents/programming-coursework/solution.txt");
    std::cout << "Solution saved to solution.txt" << std::endl;
    
    return 0;
}
