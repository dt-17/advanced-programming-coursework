#ifndef solution_visualiser_hpp
#define solution_visualiser_hpp

#include <string>
#include <vector>

class Visualiser {
public:
    // Constructor: filename of the solution file and grid dimensions.
    Visualiser(const std::string &filename, int nx, int ny, int nt);
    ~Visualiser();

    // Load the solution from file into the internal 3D vector.
    bool loadSolution();

    // Display the final time step as a 3D surface using gnuplot.
    void showFinalTimeStep();
    
    // Save the final time step image as a PNG file to the given file path.
    void saveFinalTimeStep(const std::string &outputPath);

    // Display an animation of the solution over time.
    void showAnimation(double xMin, double xMax, double yMin, double yMax, double pauseDuration);

private:
    std::string filename;
    int nx, ny, nt;

    // 3D vector to hold solution data: u[time][x][y]
    std::vector<std::vector<std::vector<double>>> u;
};

#endif // solution_visualiser_hpp
