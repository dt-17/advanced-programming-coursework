#ifndef solution_visualiser_hpp
#define solution_visualiser_hpp

#include <string>
#include <vector>

// the Visualiser class is used to load, display and save the simulation solution data
// it reads a solution file into a 3D vector and utilises GNU Plot for visualisation
class Visualiser {
public:
    // constructor: initialises the visualiser with the solution file path and grid dimensions
    // parameters:
    //   filename - path to the file containing the solution data
    //   nx       - number of grid points in the x-direction
    //   ny       - number of grid points in the y-direction
    //   nt       - number of time steps in the simulation
    Visualiser(const std::string &filename, int nx, int ny, int nt);
    
    // destructor: cleans up any resources (no explicit cleanup required as vector handles memory)
    ~Visualiser();

    // load the solution data from file into the internal 3D vector 'u'
    // returns true if the file is loaded successfully, otherwise false
    bool loadSolution();

    // display the final time step as an interactive 3D surface using GNU Plot
    // this opens a GNU Plot window to visualise the solution at the last time step
    void showFinalTimeStep();
    
    // save the final time step plot as a PNG image to the specified file path
    // parameters:
    //   outputPath - the location where the PNG image will be saved
    void saveFinalTimeStep(const std::string &outputPath);

    // display an animation of the solution evolution over time using GNU Plot
    // parameters:
    //   xMin, xMax - domain limits for the x-axis
    //   yMin, yMax - domain limits for the y-axis
    //   pauseDuration - duration (in seconds) to pause between frames
    void showAnimation(double xMin, double xMax, double yMin, double yMax, double pauseDuration, const std::string &videoOutputPath);

private:
    std::string filename;  // file path for the solution data
    int nx, ny, nt;        // grid dimensions: nx and ny for spatial, nt for temporal steps

    // 3D vector to hold the solution data, organised as u[time][x][y]
    std::vector<std::vector<std::vector<double>>> u;
};

#endif // solution_visualiser_hpp
