//
//  main.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 20/01/2025.
//

#include <iostream> // include iostream to allow inputs and outputs
#include "wave_solver.hpp"  // include solver class
#include "solution_visualiser.hpp" // include visualisation class

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

    // create an instance of the visualiser class, loading the solution data that has just been saved
    // you must pass the same filepath where you saved the solution along with the variables nx, ny and nt from above
    Visualiser vis("/Users/danieltompkins/Documents/programming-coursework/solution.txt", nx, ny, nt);
    if (!vis.loadSolution()) {
        std::cerr << "Error loading solution file." << std::endl; // return an error if data cannot be loaded
        return -1;
    }
    
    // display an interactive 3d plot of the solution at the final timestep, opens in a GNU Plot window
    //vis.showFinalTimeStep();
    
    // save a 3d plot of the solution at the final timestep in a location of your choice
    vis.saveFinalTimeStep("/Users/danieltompkins/Documents/programming-coursework/final_timestep.png");
    
    // show an animation of the solution evolution in a GNU Plot window - no outfile path is required
    // the first 4 arguments specify the domain size in x and y respectively (i.e. 0 to 10 in both)
    // final argument specifies the pause between frames in seconds (higher means greater pause and less smooth animation)
    // this can be used instead of the showFinalTimeStep() function, simply comment/uncomment the one you would like to see
    //vis.showAnimation(0.0, 10.0, 0.0, 10.0, 0.01);

    return 0;
}
